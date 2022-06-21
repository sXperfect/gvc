#!/usr/bin/env python3
import copy as cp
import logging as log
import multiprocessing as mp

import gvc.common
from . import data_structures as ds
from . import reader
from .binarization import binarize_allele_matrix, BINARIZATION_STR2ID
from .sort import sort
from .codec import CODEC_STR2ID, encode

def run_core(raw_block, ps_params, tsp_params):  
    allele_matrix, phasing_matrix, p, any_missing, not_available = raw_block
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
        any_missing,
        not_available,
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
    enc_variant = ds.GenotypePayload(new_param_set, *data_bytes)

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
    input_fpath:str,
    output_fpath,
    block_size,
    ps_params,
    tsp_params,
    num_processes
):

    lock_A = mp.Lock()
    queue_A = mp.Queue()

    lock_B = mp.Lock()
    queue_B = mp.Queue()

    procs = []

    procs.append(
        mp.Process(
            name="Reader",
            target=worker_reader, 
            args=(lock_A, queue_A, input_fpath, output_fpath, block_size),
             kwargs={'buffer_size' : num_processes*2})
    )
    
    for i_worker in range(num_processes):
        log.info('Add encoder worker {}'.format(i_worker))

        procs.append(
            mp.Process(
                name="Encoder{:02d}".format(i_worker),
                target=worker_encoder, 
                args=(lock_A, queue_A, lock_B, queue_B, ps_params, tsp_params))
        )

    procs.append(
        mp.Process(
            name="Writer",
            target=worker_writer, 
            args=(lock_B, queue_B, output_fpath), 
            kwargs={'num_processes' : num_processes})
    )

    for p in procs:
        p.start()

    for p in procs:
        p.join()

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
        dist='hrl',
        solver='nn',
        codec_name:str="jbig1",
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

        log.info('Encoding complete')


def worker_reader(lock_A, queue_A, input_fpath, output_fpath, block_size, buffer_size=1):

    if input_fpath.endswith('.txt'):
        raise RuntimeError("Deprecated!")
    elif input_fpath.endswith('.vcf') or input_fpath.endswith('vcf.gz'):
        iterator = gvc.reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
    else:
        raise ValueError('Invalid Format')

    for block_ID, raw_block in enumerate(iterator):        
        log.info('Adding block {} with size {}'.format(block_ID, raw_block[0].shape[0]))
        wait = True

        while wait:
            lock_A.acquire()

            try:
                if queue_A.qsize() < buffer_size:
                    queue_A.put(
                        [block_ID, cp.copy(raw_block)]
                    )
                    block_ID += 1

                    wait = False
            finally:
                lock_A.release()

    #? Send kill message
    wait = True
    while wait:
        lock_A.acquire()

        try:
            if queue_A.qsize() < buffer_size:
                queue_A.put(
                    [-1, None]
                )
                wait = False
        finally:
            lock_A.release()

    log.info('Stop')

def worker_encoder(lock_A, queue_A, lock_B, queue_B, ps_params, tsp_params):

    running = True
    while running:
        get_data = True
        while get_data:
            lock_A.acquire()

            try:
                if queue_A.qsize() > 0:
                    block_ID, raw_block = queue_A.get()
                    get_data = False

                    if block_ID == -1:
                        queue_A.put([-1, None])
            finally:
                lock_A.release()
        
        if block_ID == -1:

            lock_B.acquire()
            try:
                queue_B.put([-1, None, None])
            finally:
                lock_B.release()

            running = False
        else:
            block, new_param_set = run_core(raw_block, ps_params, tsp_params)            

            lock_B.acquire()
            try:
                queue_B.put(
                    [block_ID, new_param_set, block]
                )
            finally:
                lock_B.release()

    log.info('Stop')

def worker_writer(lock_B, queue_B, output_fpath, num_processes=0):
    acc_unit_param_set = None  # Act as pointer, pointing to parameter set of current AccessUnit
    acc_unit_id = 0
    blocks = []
    param_sets = []

    max_num_blocks_per_acc_unit = 2**(ds.consts.NUM_BLOCKS_LEN * 8) - 1

    curr_block_ID = 0
    block_dict = {}
    num_killed = 0

    with open(output_fpath, 'wb') as output_f:

        running = True
        retrieve_new_data = True
        while running:
        
            if retrieve_new_data:
                get_data = True
                log.info('Retrieving block')
                while get_data:
                    lock_B.acquire()

                    try:
                        if queue_B.qsize():
                            block_ID, block_param_set, block_payload = queue_B.get()

                            get_data = False

                    finally:
                        lock_B.release()

                if block_ID == -1:
                    num_killed += 1
                    log.info('Number of killed {}'.format(num_killed))

                    if num_killed >= num_processes:
                        log.info('EOF reached')

                        retrieve_new_data = False
                        # break
                else:
                    log.info('Storing block')
                    block_dict[block_ID] = [block_param_set, block_payload]

            if curr_block_ID in block_dict:
                log.info('Writing block {}'.format(curr_block_ID))

                new_param_set, block = block_dict[curr_block_ID]
                del block_dict[curr_block_ID]
                curr_block_ID += 1
                
                # If parameter set of current block different from parameter set of current access unit,
                # store blocks as access unit
                if acc_unit_param_set is None:
                    log.info('Parameter set of access unit is None -> set to new parameter set')
                    acc_unit_param_set = new_param_set
                    param_sets.append(acc_unit_param_set)
                    output_f.write(acc_unit_param_set.to_bytes())

                elif new_param_set != acc_unit_param_set or len(blocks) == max_num_blocks_per_acc_unit:
                    log.info('Parameter set and parameter set of current access unit is different')

                    # Store blocks as an Access Unit
                    log.info('Store access unit ID {:03d}, num blocks: {:d}'.format(acc_unit_id, len(blocks)))
                    gvc.common.store_access_unit(output_f, acc_unit_id, acc_unit_param_set, blocks)

                    # Initialize values for the new AccessUnit
                    acc_unit_id += 1
                    blocks.clear()

                    # Check if similar parameter set is already created before                       
                    is_param_set_unique = True
                    for stored_param_set in param_sets:
                        if stored_param_set == new_param_set:
                            is_param_set_unique = False
                            break
                    
                    # If parameter set is unique, store in list of parameter sets and store in GVC file
                    if is_param_set_unique:
                        log.info('New parameter set is unique')
                        new_param_set.parameter_set_id = len(param_sets)

                        acc_unit_param_set = new_param_set

                        param_sets.append(acc_unit_param_set)
                        output_f.write(acc_unit_param_set.to_bytes())

                    else:
                        log.info('New parameter set is not unique')
                        del new_param_set
                        acc_unit_param_set = stored_param_set

                blocks.append(block)

            if len(block_dict) == 0 and not retrieve_new_data:
                log.info("No new block available")
                running = False

        if len(blocks):
            lock_B.acquire()

            try:
                qsize = queue_B.qsize()

            finally:
                lock_B.release()

            log.info("QueueB qsize: {}".format(qsize))
            log.info("Dict: {}".format(block_dict))
            log.info("Dict size: {}".format(len(block_dict)))
            log.info("curr_block_ID:{}".format(curr_block_ID))

            # Store the remaining blocks
            log.info('Store the remaining blocks')
            log.info('Store access unit ID {:03d}, num blocks: {:d}'.format(acc_unit_id, len(blocks)))
            gvc.common.store_access_unit(output_f, acc_unit_id, acc_unit_param_set, blocks)

    log.info('Stop')