#!/usr/bin/env python3
import typing as t
import copy as cp
import logging as log
import multiprocessing as mp

import gvc.common
from . import data_structures as ds
from . import reader
from .binarization import binarize_allele_matrix, BINARIZATION_STR2ID
from .sort import sort
from .codec import CODEC_STR2ID, encode

def run_core(
    raw_block:t.List, 
    ps_params:t.List,
    tsp_params:t.List
):
    
    allele_matrix, phasing_matrix, p, missing_rep_val, na_rep_val = raw_block
    binarization_id,codec_id,axis,sort_rows,sort_cols,transpose = ps_params
    dist_f_name, solver_name, solver_profile = tsp_params

    # Execute part 4.2 - binarization of allele matrix
    log.info('Execute part 4.2 - Binarization')
    bin_allele_matrices, additional_info = binarize_allele_matrix(
        allele_matrix, 
        binarization_id, 
        axis=axis
    )
                
    #? Create parameter based on binarization and encoder parameter
    log.info('Create parameter set')
    new_param_set = gvc.common.create_parameter_set(
        missing_rep_val,
        na_rep_val,
        p,
        phasing_matrix,
        additional_info,
        binarization_id,
        codec_id,
        axis,
        sort_rows,
        sort_cols,
        transpose=transpose
    )

    # Execute part 4.3 - sorting
    log.info('Execute part 4.3 - Sorting')
    sorted_data = sort(
        new_param_set, 
        bin_allele_matrices, 
        phasing_matrix, 
        dist_f_name=dist_f_name, 
        solver_name=solver_name,
        solver_profile=solver_profile,
    )

    # Execute part 4.4 - entropy coding
    log.info('Execute part 4.4')
    data_bytes = encode(new_param_set, additional_info, *sorted_data)

    # Initialize EncodedVariant, ParameterSet is not stored internally in EncodedVariants
    # Order of arguments here is important (See EncodedVariants)
    enc_variant = ds.GenotypePayload(new_param_set, *data_bytes, missing_rep_val, na_rep_val)

    # Create new Block and store
    block = ds.Block.from_encoded_variant(enc_variant)
    
    return block, new_param_set

def run_no_threads(
    input_fpath:str,
    output_fpath,
    block_size,
    ps_params,
    tsp_params,
):
    log.info('run without multithreading')

    with open(output_fpath, 'wb') as output_f:

        ac_unit_param_set = None  # Act as pointer, pointing to parameter set of current AccessUnit
        acc_unit_id = 0
        blocks = []
        param_sets = []

        num_bytes_per_block = []
        max_num_blocks_per_acc_unit = 2**(ds.consts.NUM_BLOCKS_LEN * 8) - 1

        if input_fpath.endswith('.vcf'):
            iterator = reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
        elif input_fpath.endswith('vcf.gz'):
            iterator = reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
        else:
            raise ValueError('Invalid Format')

        for block_ID, raw_block in enumerate(iterator):
            log.info(f"Processing block {block_ID}")          
            block, new_param_set = run_core(raw_block, ps_params, tsp_params)
            
            num_bytes_per_block.append(len(block))
            
            #? If parameter set of current block different from parameter set of current access unit,
            #? store blocks as access unit
            if ac_unit_param_set is None:
                log.info('Set to new parameter set')
                ac_unit_param_set = new_param_set
                param_sets.append(ac_unit_param_set)
                output_f.write(ac_unit_param_set.to_bytes())

            elif new_param_set != ac_unit_param_set or len(blocks) == max_num_blocks_per_acc_unit:
                log.info('New parameter set found, store blocks')

                # Store blocks as an Access Unit
                log.info('Store access unit ID {:03d}'.format(acc_unit_id))
                gvc.common.store_access_unit(output_f, acc_unit_id, ac_unit_param_set, blocks)

                # Initialize values for the new AccessUnit
                acc_unit_id += 1
                blocks.clear()

                #? Check if similar parameter set is already created before                       
                is_param_set_unique = True
                for stored_param_set in param_sets:
                    if stored_param_set == new_param_set:
                        is_param_set_unique = False
                        break
                
                #? If parameter set is unique, store in list of parameter sets and store in GVC file
                if is_param_set_unique:
                    log.info('New parameter set is unique')
                    new_param_set.parameter_set_id = len(param_sets)

                    ac_unit_param_set = new_param_set

                    param_sets.append(ac_unit_param_set)
                    output_f.write(ac_unit_param_set.to_bytes())

                else:
                    log.info('New parameter set is not unique')
                    del new_param_set
                    ac_unit_param_set = stored_param_set

            blocks.append(block)

        if len(blocks):
            #? Store the remaining blocks
            log.info('Store the remaining blocks')
            gvc.common.store_access_unit(output_f, acc_unit_id, ac_unit_param_set, blocks)
            
def run_multiprocessing(
    input_fpath,
    output_fpath,
    block_size,
    ps_params,
    tsp_params,
    num_processes,
):
    if not isinstance(num_processes, int) or isinstance(num_processes, bool):
        raise TypeError("num_processes must be an integer")
    if num_processes < 1:
        raise ValueError("num_processes must be positive")

    queue_a = mp.Queue(maxsize=max(1, num_processes * 2))
    queue_b = mp.Queue(maxsize=max(1, num_processes * 2))

    reader_proc = mp.Process(
        name="Reader",
        target=worker_reader,
        args=(queue_a, input_fpath, output_fpath, block_size, num_processes),
    )
    encoder_procs = [
        mp.Process(
            name="Encoder{:02d}".format(i_worker),
            target=worker_encoder,
            args=(queue_a, queue_b, ps_params, tsp_params),
        )
        for i_worker in range(num_processes)
    ]
    writer_proc = mp.Process(
        name="Writer",
        target=worker_writer,
        args=(queue_b, output_fpath, num_processes),
    )

    procs = [reader_proc] + encoder_procs + [writer_proc]
    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()

    failures = [
        "{} exited with status {}".format(proc.name, proc.exitcode)
        for proc in procs
        if proc.exitcode != 0
    ]
    queue_a.close()
    queue_b.close()
    if failures:
        raise RuntimeError(
            "multiprocessing encoder failed: {}".format("; ".join(failures))
        )

class Encoder(object):

    def __init__(self,
        input_fpath,
        output_fpath,
        binarization_name:str="bit_plane",
        axis:int=1,
        sort_cols=True,
        sort_rows=True,
        transpose=False,
        block_size=2536,
        max_cols=None,
        dist='ham',
        solver='nn',
        codec_name:str="jbig",
        preset_mode=1,
        num_threads=0,
    ):

        self.input_fpath = input_fpath
        self.output_fpath = output_fpath

        # Parameter Set
        self.binarization_id = BINARIZATION_STR2ID[binarization_name]
        self.codec_id = CODEC_STR2ID[codec_name]
        self.axis = axis
        self.sort_cols = sort_cols
        self.sort_rows = sort_rows
        self.transpose = transpose

        # Binarization parameter (additional)
        self.block_size = block_size
        self.max_cols = max_cols
        
        # Parameter for sorting process
        self.dist = dist
        self.solver = solver

        # Additional parameter
        self.preset_mode = preset_mode
        self.num_threads = num_threads
        
    @property
    def ps_params(self):
        return [
            self.binarization_id,
            self.codec_id,
            self.axis,
            self.sort_rows,
            self.sort_cols,
            self.transpose,
        ]
        
    @property
    def tsp_params(self):
        return [
            self.dist,
            self.solver,
            self.preset_mode
        ]

    def run(self):
        log.info('encoding: {} -> {}'.format(self.input_fpath, self.output_fpath))

        if self.num_threads == 0:
            run_no_threads(
                self.input_fpath,
                self.output_fpath,
                self.block_size,
                self.ps_params,
                self.tsp_params,
            )

        elif self.num_threads > 0:
            run_multiprocessing(
                self.input_fpath,
                self.output_fpath,
                self.block_size,
                self.ps_params,
                self.tsp_params,
                self.num_threads,
            )

        else:
            log.error('Invalid value for num_threads')
            raise ValueError('Invalid value for num_threads')

        log.debug('Encoding complete')


def worker_reader(
    queue_a,
    input_fpath,
    output_fpath,
    block_size,
    num_processes,
):
    if not (
        input_fpath.endswith(".vcf")
        or input_fpath.endswith(".vcf.gz")
    ):
        raise ValueError("Invalid Format")

    try:
        iterator = gvc.reader.vcf_genotypes_reader(
            input_fpath, output_fpath, block_size
        )
        for block_id, raw_block in enumerate(iterator):
            log.info(
                "Adding block {} with size {}".format(
                    block_id, raw_block[0].shape[0]
                )
            )
            queue_a.put((block_id, cp.copy(raw_block)))
    finally:
        # One sentinel per worker: no worker has to put a sentinel back into
        # the queue, which avoids races and makes shutdown deterministic.
        for _ in range(num_processes):
            queue_a.put(None)


def worker_encoder(queue_a, queue_b, ps_params, tsp_params):
    while True:
        item = queue_a.get()
        if item is None:
            queue_b.put(None)
            return

        block_id, raw_block = item
        block, new_param_set = run_core(raw_block, ps_params, tsp_params)
        queue_b.put((block_id, new_param_set, block))


def _store_ordered_block(
    output_f,
    block,
    new_param_set,
    state,
):
    max_num_blocks = 2 ** (ds.consts.NUM_BLOCKS_LEN * 8) - 1
    ac_unit_param_set = state["parameter_set"]
    blocks = state["blocks"]
    param_sets = state["parameter_sets"]

    if ac_unit_param_set is None:
        ac_unit_param_set = new_param_set
        param_sets.append(ac_unit_param_set)
        output_f.write(ac_unit_param_set.to_bytes())

    elif new_param_set != ac_unit_param_set or len(blocks) == max_num_blocks:
        gvc.common.store_access_unit(
            output_f,
            state["access_unit_id"],
            ac_unit_param_set,
            blocks,
        )
        state["access_unit_id"] += 1
        blocks.clear()

        stored_match = next(
            (
                stored
                for stored in param_sets
                if stored == new_param_set
            ),
            None,
        )
        if stored_match is None:
            new_param_set.parameter_set_id = len(param_sets)
            ac_unit_param_set = new_param_set
            param_sets.append(ac_unit_param_set)
            output_f.write(ac_unit_param_set.to_bytes())
        else:
            ac_unit_param_set = stored_match

    blocks.append(block)
    state["parameter_set"] = ac_unit_param_set


def worker_writer(queue_b, output_fpath, num_processes):
    if num_processes < 1:
        raise ValueError("num_processes must be positive")

    state = {
        "parameter_set": None,
        "access_unit_id": 0,
        "blocks": [],
        "parameter_sets": [],
    }
    pending = {}
    next_block_id = 0
    stopped_workers = 0

    with open(output_fpath, "wb") as output_f:
        while stopped_workers < num_processes:
            item = queue_b.get()
            if item is None:
                stopped_workers += 1
                continue

            block_id, new_param_set, block = item
            if block_id in pending or block_id < next_block_id:
                raise ValueError(
                    "duplicate or stale encoded block id {}".format(block_id)
                )
            pending[block_id] = (new_param_set, block)

            while next_block_id in pending:
                current_param_set, current_block = pending.pop(next_block_id)
                _store_ordered_block(
                    output_f,
                    current_block,
                    current_param_set,
                    state,
                )
                next_block_id += 1

        if pending:
            raise RuntimeError(
                "missing encoded block before block {}".format(
                    min(pending)
                )
            )

        if state["blocks"]:
            gvc.common.store_access_unit(
                output_f,
                state["access_unit_id"],
                state["parameter_set"],
                state["blocks"],
            )

    log.info("Stop")

