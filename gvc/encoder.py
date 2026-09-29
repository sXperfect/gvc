#!/usr/bin/env python3
import typing as t
import copy as cp
import logging as log
import multiprocessing as mp
import os
import pickle
import queue
import uuid

import gvc.common
from . import data_structures as ds
from . import reader
from .binarization import binarize_allele_matrix, BINARIZATION_STR2ID
from .sort import sort
from .codec import CODEC_STR2ID, encode
from .multiprocessing import (
    EncodedBlock,
    EncodeProcessSupervisor,
    Progress,
    ReaderDone,
    StopWork,
    WorkItem,
    WorkerDone,
    WriterDone,
)
from .multiprocessing.supervisor import child_entry

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
            
def _queue_put(target_queue, item, stop_event, timeout=0.2):
    while not stop_event.is_set():
        try:
            target_queue.put(item, timeout=timeout)
            return True
        except queue.Full:
            continue
    return False


def _queue_get(source_queue, stop_event, timeout=0.2):
    while not stop_event.is_set():
        try:
            return source_queue.get(timeout=timeout)
        except queue.Empty:
            continue
    return None


def _validate_parallel_output_path(output_fpath):
    output_path = os.path.abspath(output_fpath)
    parent = os.path.dirname(output_path) or os.getcwd()
    if not os.path.exists(parent):
        raise FileNotFoundError(
            "output directory does not exist: {}".format(parent)
        )
    if not os.path.isdir(parent):
        raise NotADirectoryError(
            "output parent is not a directory: {}".format(parent)
        )
    if os.path.isdir(output_path):
        raise IsADirectoryError(
            "output path is a directory: {}".format(output_path)
        )

    metadata_path = output_path + ".metadata"
    if os.path.exists(metadata_path) and not os.path.isdir(metadata_path):
        raise NotADirectoryError(
            "metadata sidecar path is not a directory: {}".format(
                metadata_path
            )
        )
    return output_path


def _temp_output_path(output_fpath):
    return "{}.tmp.{}.{}".format(
        output_fpath,
        os.getpid(),
        uuid.uuid4().hex,
    )


class _BlockProcessingError(RuntimeError):
    def __init__(self, block_id, message):
        self.block_id = block_id
        super().__init__(message)


def run_multiprocessing(
    input_fpath,
    output_fpath,
    block_size,
    ps_params,
    tsp_params,
    num_processes,
    start_method=None,
    stall_timeout=None,
    process_initializer=None,
    process_initializer_args=(),
):
    if not isinstance(num_processes, int) or isinstance(num_processes, bool):
        raise TypeError("num_processes must be an integer")
    if num_processes < 1:
        raise ValueError("num_processes must be positive")

    output_fpath = _validate_parallel_output_path(output_fpath)

    if process_initializer is not None and not callable(process_initializer):
        raise TypeError("process_initializer must be callable or None")
    if process_initializer_args is None:
        process_initializer_args = ()
    else:
        process_initializer_args = tuple(process_initializer_args)

    if start_method == "spawn" and process_initializer is not None:
        try:
            pickle.dumps((process_initializer, process_initializer_args))
        except Exception as exc:
            raise TypeError(
                "spawn process_initializer and arguments must be picklable"
            ) from exc

    context = (
        mp.get_context(start_method)
        if start_method is not None
        else mp.get_context()
    )
    work_q = context.Queue(maxsize=max(1, num_processes * 2))
    result_q = context.Queue(maxsize=max(1, num_processes * 2))
    error_q = context.Queue()
    status_q = context.Queue()
    stop_event = context.Event()

    temp_output = _temp_output_path(output_fpath)

    reader_proc = context.Process(
        name="GVC-Reader",
        target=child_entry,
        args=(
            "reader",
            None,
            error_q,
            stop_event,
            worker_reader,
            (
                work_q,
                status_q,
                stop_event,
                input_fpath,
                temp_output,
                block_size,
                num_processes,
            ),
            None,
            (),
        ),
    )
    reader_proc._gvc_stage = "reader"
    reader_proc._gvc_worker_id = None

    encoder_procs = []
    for worker_id in range(num_processes):
        proc = context.Process(
            name="GVC-Encoder{:02d}".format(worker_id),
            target=child_entry,
            args=(
                "encoder",
                worker_id,
                error_q,
                stop_event,
                worker_encoder,
                (
                    worker_id,
                    work_q,
                    result_q,
                    status_q,
                    stop_event,
                    ps_params,
                    tsp_params,
                ),
                process_initializer,
                process_initializer_args,
            ),
        )
        proc._gvc_stage = "encoder"
        proc._gvc_worker_id = worker_id
        encoder_procs.append(proc)

    writer_proc = context.Process(
        name="GVC-Writer",
        target=child_entry,
        args=(
            "writer",
            None,
            error_q,
            stop_event,
            worker_writer,
            (
                result_q,
                status_q,
                stop_event,
                temp_output,
                num_processes,
            ),
            None,
            (),
        ),
    )
    writer_proc._gvc_stage = "writer"
    writer_proc._gvc_worker_id = None

    supervisor = EncodeProcessSupervisor(
        processes=[reader_proc] + encoder_procs + [writer_proc],
        error_q=error_q,
        status_q=status_q,
        stop_event=stop_event,
        queues=[work_q, result_q],
        temp_output=temp_output,
        final_output=output_fpath,
        stall_timeout=stall_timeout,
    )
    return supervisor.run()

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
        multiprocessing_start_method=None,
        multiprocessing_stall_timeout=None,
        multiprocessing_initializer=None,
        multiprocessing_initializer_args=(),
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
        self.multiprocessing_start_method = multiprocessing_start_method
        self.multiprocessing_stall_timeout = multiprocessing_stall_timeout
        self.multiprocessing_initializer = multiprocessing_initializer
        self.multiprocessing_initializer_args = tuple(
            multiprocessing_initializer_args or ()
        )
        
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
                start_method=self.multiprocessing_start_method,
                stall_timeout=self.multiprocessing_stall_timeout,
                process_initializer=self.multiprocessing_initializer,
                process_initializer_args=self.multiprocessing_initializer_args,
            )

        else:
            log.error('Invalid value for num_threads')
            raise ValueError('Invalid value for num_threads')

        log.debug('Encoding complete')


def worker_reader(
    work_q,
    status_q,
    stop_event,
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

    iterator = gvc.reader.vcf_genotypes_reader(
        input_fpath, output_fpath, block_size
    )
    total_blocks = 0
    for block_id, raw_block in enumerate(iterator):
        if stop_event.is_set():
            return
        item = WorkItem(block_id, cp.copy(raw_block))
        if not _queue_put(work_q, item, stop_event):
            return
        total_blocks += 1
        _queue_put(
            status_q,
            Progress("reader", total_blocks),
            stop_event,
        )

    for _ in range(num_processes):
        if not _queue_put(work_q, StopWork(), stop_event):
            return
    _queue_put(status_q, ReaderDone(total_blocks), stop_event)


def worker_encoder(
    worker_id,
    work_q,
    result_q,
    status_q,
    stop_event,
    ps_params,
    tsp_params,
):
    processed_blocks = 0
    while not stop_event.is_set():
        item = _queue_get(work_q, stop_event)
        if item is None:
            return
        if isinstance(item, StopWork):
            _queue_put(result_q, WorkerDone(worker_id), stop_event)
            return
        if not isinstance(item, WorkItem):
            raise TypeError(
                "unexpected work-queue message: {}".format(type(item).__name__)
            )

        try:
            block, new_param_set = run_core(
                item.raw_block, ps_params, tsp_params
            )
        except BaseException as exc:
            raise _BlockProcessingError(
                item.block_id, str(exc)
            ) from exc

        if not _queue_put(
            result_q,
            EncodedBlock(item.block_id, new_param_set, block),
            stop_event,
        ):
            return
        processed_blocks += 1
        _queue_put(
            status_q,
            Progress(
                "encoder[{}]".format(worker_id),
                processed_blocks,
            ),
            stop_event,
        )


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
            (stored for stored in param_sets if stored == new_param_set),
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


def worker_writer(
    result_q,
    status_q,
    stop_event,
    output_fpath,
    num_processes,
):
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
    written_blocks = 0

    with open(output_fpath, "wb") as output_f:
        while (
            stopped_workers < num_processes
            and not stop_event.is_set()
        ):
            item = _queue_get(result_q, stop_event)
            if item is None:
                return

            if isinstance(item, WorkerDone):
                stopped_workers += 1
                continue
            if not isinstance(item, EncodedBlock):
                raise TypeError(
                    "unexpected result-queue message: {}".format(
                        type(item).__name__
                    )
                )

            block_id = item.block_id
            if block_id in pending or block_id < next_block_id:
                raise ValueError(
                    "duplicate or stale encoded block id {}".format(block_id)
                )
            pending[block_id] = (item.parameter_set, item.block)

            while next_block_id in pending:
                current_param_set, current_block = pending.pop(next_block_id)
                _store_ordered_block(
                    output_f,
                    current_block,
                    current_param_set,
                    state,
                )
                next_block_id += 1
                written_blocks += 1
                _queue_put(
                    status_q,
                    Progress("writer", written_blocks),
                    stop_event,
                )

        if stop_event.is_set():
            return
        if pending:
            raise RuntimeError(
                "missing encoded block before block {}".format(min(pending))
            )

        if state["blocks"]:
            gvc.common.store_access_unit(
                output_f,
                state["access_unit_id"],
                state["parameter_set"],
                state["blocks"],
            )
        output_f.flush()
        os.fsync(output_f.fileno())

    _queue_put(status_q, WriterDone(written_blocks), stop_event)
    log.info("Stop")

