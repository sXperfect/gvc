#!/usr/bin/env python3

import logging as log

import gvc.reader
from gvc.reader import FORMAT_ID
import gvc.common
import gvc.processes
from . import data_structures as ds
import gvc.binarization
from . import reader
from .sort import sort
from .codec import encode

def _run_no_threads(
    # input_f,
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
    solver_profile
):

    log.info('run without multithreading')

    with open(output_fpath, 'wb') as output_f:

        ac_unit_param_set = None  # Act as pointer, pointing to parameter set of current AccessUnit
        acc_unit_id = 0
        blocks = []
        param_sets = []

        num_bytes_per_block = []

        max_num_blocks_per_acc_unit = 2**(ds.consts.NUM_BLOCKS_LEN * 8) - 1

        lines = []

        if input_fpath.endswith('.txt'):
            raise RuntimeError("Deprecated!")
        elif input_fpath.endswith('.vcf'):
            iterator = gvc.reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
        elif input_fpath.endswith('vcf.gz'):
            iterator = gvc.reader.vcf_genotypes_reader(input_fpath, output_fpath, block_size)
        else:
            raise ValueError('Invalid Format')

        for block_ID, raw_block in enumerate(iterator):
            allele_matrix, phasing_matrix, p, any_missing, not_available = raw_block
            
            log.info('Adding block {} with size {}'.format(block_ID, allele_matrix.shape[0]))
            # block, num_samples, p = raw_block
            # format_id = raw_block[0]
            
            # if format_id == FORMAT_ID.VCF:
            #     allele_matrix, phasing_matrix, p, any_missing, not_available = reader.raw_vcf_to_gt_mat(raw_block[1:])
            # else:
            #     raise NotImplementedError("Not yet implemented or deprecated!")

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

            num_bytes_per_block.append(len(block))
            
            # If parameter set of current block different from parameter set of current access unit,
            # store blocks as access unit
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

                    ac_unit_param_set = new_param_set

                    param_sets.append(ac_unit_param_set)
                    output_f.write(ac_unit_param_set.to_bytes())

                else:
                    log.info('New parameter set is not unique')
                    del new_param_set
                    ac_unit_param_set = stored_param_set

            blocks.append(block)

            lines.clear()

        if len(blocks):
            # Store the remaining blocks
            log.info('Store the remaining blocks')
            gvc.common.store_access_unit(output_f, acc_unit_id, ac_unit_param_set, blocks)

class Encoder(object):

    def __init__(self,
        input_fpath,
        output_fpath,
        binarization="bit_plane",
        axis=1,
        sort_cols=True,
        sort_rows=True,
        transpose=False,
        block_size=256,
        max_cols=None,
        dist='hrl',
        solver='nn',
        encoder="jbig1",
        preset_mode=1,
        num_threads=0,
    ):

        self.input_fpath = input_fpath
        self.output_fpath = output_fpath

        # Parameter Set
        self.binarization = binarization
        self.encoder = encoder
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

    def run(self):
        log.info('encoding: {} -> {}'.format(self.input_fpath, self.output_fpath))

        if self.num_threads == 0:
            _run_no_threads(
                self.input_fpath,
                self.output_fpath,
                self.binarization,
                self.encoder,
                self.axis,
                self.sort_rows,
                self.sort_cols,
                self.transpose,
                self.block_size,
                self.dist,
                self.solver,
                self.preset_mode
            )

        elif self.num_threads > 0:
            # input_f = open(self.input_fpath, 'r')
            gvc.processes.run_multiprocessing(
                self.input_fpath,
                self.output_fpath,
                self.binarization,
                self.encoder,
                self.axis,
                self.sort_rows,
                self.sort_cols,
                self.transpose,
                self.block_size,
                self.dist,
                self.solver,
                self.preset_mode,
                self.num_threads
            )

        else:
            log.error('Invalid value for num_threads')
            raise ValueError('Invalid value for num_threads')

        log.info('Encoding complete')


