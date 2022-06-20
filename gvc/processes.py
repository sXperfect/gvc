import copy as cp
import logging as log
import multiprocessing as mp

import numpy as np

import gvc.common

from . import data_structures as ds
from .reader import FORMAT_ID, raw_vcf_genotypes_reader, raw_vcf_to_gt_mat
from .sort import sort
from .codec import encode

def worker_reader(lock_A, queue_A, input_fpath, output_fpath, block_size, buffer_size=1):

    if input_fpath.endswith('.txt'):
        raise RuntimeError("Deprecated!")
    elif input_fpath.endswith('.vcf') or input_fpath.endswith('vcf.gz'):
        iterator = gvc.reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
    else:
        raise ValueError('Invalid Format')

    for block_ID, raw_block in enumerate(iterator):
        allele_matrix, phasing_matrix, p, any_missing, not_available = raw_block
        
        log.info('Adding block {} with size {}'.format(block_ID, allele_matrix.shape[0]))

        # format_id = raw_block[0]
        
        # if format_id == FORMAT_ID.TXT:
        #     raise RuntimeError("Deprecated!")
        #     allele_matrix = raw_block[0]
        #     log.info('Adding block {} with size {}'.format(block_ID, allele_matrix.shape[0]))
        # elif format_id == FORMAT_ID.VCF:
        #     log.info('Adding block {} with size {}'.format(block_ID, len(raw_block[1])))

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

def worker_encoder(lock_A, queue_A, lock_B, queue_B, param_set_list):

    binarizer, encoder, axis, sort_rows, sort_cols, transpose, dist_f_name, solver_name, solver_profile = param_set_list

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
            # format_id = raw_block[0]
            
            # if format_id == FORMAT_ID.VCF:
            #     allele_matrix, phasing_matrix, p, any_missing, not_available = raw_vcf_to_gt_mat(raw_block[1:])
            # else:
            #     raise NotImplementedError("Not yet implemented or deprecated!")
            
            allele_matrix, phasing_matrix, p, any_missing, not_available = raw_block

            log.info('Encoding block {} with size {}'.format(block_ID, allele_matrix.shape[0]))

            # Execute part 4.2 - binarization of allele matrix
            log.info('Execute part 4.2 - Binarization')
            bin_allele_matrices, additional_info = gvc.binarization.binarize_allele_matrix(
                allele_matrix, 
                binarizer=binarizer, 
                axis=axis
            )
            
            # Create parameter based on binarization and encoder parameter
            log.info('Create parameter set')
            new_param_set = gvc.common.create_parameter_set(
                any_missing,
                not_available,
                p,
                phasing_matrix,
                additional_info,
                binarizer,
                encoder,
                axis,
                sort_rows,
                sort_cols,
                transpose,
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


def run_multiprocessing(
    input_fpath:str,
    output_fpath,
    binarizer,
    encoder,
    axis,
    sort_rows,
    sort_cols,
    transpose,
    block_size,
    dist_f_name,
    solver_name,
    solver_profile,
    num_processes
):

    lock_A = mp.Lock()
    queue_A = mp.Queue()

    lock_B = mp.Lock()
    queue_B = mp.Queue()

    param_set_list = [binarizer, encoder, axis, sort_rows, sort_cols, transpose, dist_f_name, solver_name, solver_profile]

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
                args=(lock_A, queue_A, lock_B, queue_B, param_set_list))
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

