import time
import logging as log

import numpy as np
# from ds.data_unit import DataUnitHeader

import gvc.utils
import gvc.data_structures
import gvc.bitstream
import gvc.binarization
from .data_structures import consts
from . import data_structures as ds

from . import codec
from .codec import jbig
from . import cdebinarize
from . import cquery

def _decode_encoded_variants(
    param_set:ds.ParameterSet,
    encoded_variants, 
):

    bin_matrices = np.empty(param_set.num_variants_flags, dtype=object)
    for i_bin_mat in range(param_set.num_variants_flags):

        bin_mat = codec.decode(
            encoded_variants.variants_payloads[i_bin_mat],
            encoded_variants.variants_row_ids_payloads[i_bin_mat],
            encoded_variants.variants_col_ids_payloads[i_bin_mat],
            param_set.variants_coder_ids[i_bin_mat],
        )
        
        if param_set.transpose_variants_mat_flags[i_bin_mat]:
            log.info('Re-transpose binary matrix')
            bin_mat = bin_mat.T

        bin_matrices[i_bin_mat] = bin_mat

    if param_set.binarization_flag == consts.BinarizationID.BIT_PLANE:
        additional_info = param_set.num_bin_mat
    elif param_set.binarization_flag in [consts.BinarizationID.ROW_SPLIT, consts.BinarizationID.ROW_BIN_SPLIT]:
        additional_info = ds.VectorAMax.from_bytes(encoded_variants.variants_amax_payload.read()).vector
    elif param_set.binarization_flag == consts.BinarizationID.ROW_BIN_SPLIT2:
        additional_info = None
    else:
        raise ValueError()

    allele_matrix = gvc.binarization.debinarize_bin_matrices(
        bin_matrices, additional_info, param_set.concat_axis, param_set.binarization_flag
    )

    allele_matrix = gvc.binarization.undo_adaptive_max_value(
        allele_matrix, param_set.any_missing_flag, param_set.not_available_flag)

    if param_set.encode_phase_data:
        phasing_matrix = codec.decode(
            encoded_variants.phase_payload,
            encoded_variants.phase_row_ids_payload,
            encoded_variants.phase_col_ids_payload,
            param_set.phase_coder_ids,
        )
        
        if param_set.transpose_phase_mat_flag:
            log.info('Re-transpose phase matrix')
            phasing_matrix = phasing_matrix.T
    else:
        phasing_matrix = param_set.phase_value

    return gvc.binarization.reconstruct_genotype_matrix(allele_matrix, phasing_matrix, param_set.p)
    # return allele_matrix, phasing_matrix

def _compute_amax(amax_vec):
    cumsum_amax_vec = np.cumsum(amax_vec)

    return cumsum_amax_vec, amax_vec


def _selective_decode_encoded_variants(
    param_set:ds.ParameterSet, 
    encoded_variants, 
    row_slice, 
    query_col_ids,
):
    # TODO: Column mask

    if param_set.binarization_flag == gvc.binarization.BIT_PLANE_FLAG:
        bin_matrices = np.empty(param_set.num_variants_flags, dtype=object)
        for i_bin_mat in range(param_set.num_variants_flags):

            bin_mat = codec.decode(
                encoded_variants.variants_payloads[i_bin_mat],
                encoded_variants.variants_row_ids_payloads[i_bin_mat],
                encoded_variants.variants_col_ids_payloads[i_bin_mat],
                param_set.variants_coder_ids[i_bin_mat],
            )
            
            if query_col_ids is None:
                if col_ids is not None:
                    bin_mat = bin_mat[:, col_ids]
                else:
                    pass
            else:
                nqci = cquery.cget_col_ids(query_col_ids, param_set.p)
                
                try:
                    col_ids = col_ids[nqci]
                except:
                    col_ids = nqci
                    
                bin_mat = bin_mat[:, col_ids]  
            
            if param_set.transpose_variants_mat_flags[i_bin_mat]:
                log.info('Re-transpose binary matrix')
                bin_mat = bin_mat.T

            bin_matrices.append(bin_mat)
            
        raise NotImplementedError()
        additional_info = param_set.num_bin_mat
    elif param_set.binarization_flag in [consts.BinarizationID.ROW_BIN_SPLIT]:       
        bin_mat, row_ids, col_ids = codec.decode(
            encoded_variants.variants_payloads[0],
            encoded_variants.variants_row_ids_payloads[0],
            encoded_variants.variants_col_ids_payloads[0],
            param_set.variants_coder_ids[0],
            unsort=False
        )
        
        if query_col_ids is None:
            if col_ids is not None:
                bin_mat = bin_mat[:, col_ids]
            else:
                pass
        else:
            nqci = cquery.cget_col_ids(query_col_ids, param_set.p)
            
            try:
                col_ids = col_ids[nqci]
            except:
                col_ids = nqci
                
            bin_mat = bin_mat[:, col_ids]            
            
        amax_vec = ds.VectorAMax.from_bytes(
            encoded_variants.variants_amax_payload.read()
        ).vector

        selected_row_ids = row_ids[row_slice]

        min_row_ids = np.min(selected_row_ids)
        max_row_ids = np.max(selected_row_ids)

        bin_mat_slice = slice(
            np.sum(amax_vec[:min_row_ids]),
            np.sum(amax_vec[:max_row_ids+1]),
        )

        if param_set.binarization_flag ==  consts.BinarizationID.ROW_BIN_SPLIT:
            allele_matrix = gvc.binarization.debin_row_bin_split(bin_mat[bin_mat_slice, :], amax_vec[min_row_ids:max_row_ids+1])
            selected_allele_matrix = allele_matrix[selected_row_ids-min_row_ids, :]
        else:
            raise NotImplementedError()
    else:
        raise ValueError("Invalid binarization flag")

    if param_set.encode_phase_data:
        phasing_matrix = codec.decode(
            encoded_variants.phase_payload,
            encoded_variants.phase_row_ids_payload,
            encoded_variants.phase_col_ids_payload,
            param_set.phase_coder_ids,
        )

        # TODO
        raise NotImplementedError("Use cheaper method for selective decoding")
        
        if param_set.transpose_phase_mat_flag:
            log.info('Re-transpose phase matrix')
            phasing_matrix = phasing_matrix.T
    else:
        phase_val = param_set.phase_value

        gt_mat = cdebinarize.reconstruct_genotype_matrix_using_phase_value(
            selected_allele_matrix, phase_val, param_set.p
        )

        return gt_mat

def _decode_block_payload(param_set, block: ds.Block):
    log.info('decoding block payload')

    if block.block_header.content_id == consts.ContentID.GENOTYPE:
        log.info('Decoding EncodedVariants')

        return _decode_encoded_variants(param_set, block.block_payload)

    else:
        error_msg = 'Invalid content id: {}'.format(block.block_header.content_id)
        log.error(error_msg)
        raise gvc.errors.GvcError(error_msg)


def _decode_access_unit(decoder_context):
    log.info('Decoding access unit')

    allele_tensors = []
    phasing_tensors = []
    for i, block in enumerate(decoder_context.curr_access_unit.blocks):
        log.info('Decoding block {}'.format(i))
        allele_tensor, phasing_tensor = _decode_block_payload(decoder_context.curr_parameter_set, block)

        allele_tensors.append(allele_tensor)
        phasing_tensors.append(phasing_tensor)

    return np.concatenate(allele_tensors, axis=0), np.concatenate(phasing_tensors, axis=0)

def _decode_and_write_access_unit(out_f, decoder_context):
    log.info('Decoding access unit')

    for i_block, block in enumerate(decoder_context.curr_access_unit.blocks):
        log.info('Decoding block {}'.format(i_block))

        stime = time.time()
        allele_tensor, phasing_tensor = _decode_block_payload(decoder_context.curr_parameter_set, block)

        lines = gvc.binarization.tensor_to_txt(allele_tensor, phasing_tensor)
        nrows, ncols = lines.shape
        for i in range(nrows):
            out_f.write("\t".join(lines[i, :]) + "\n")

            for j in range(ncols):
                out_f.write(lines[i,j])

            if j < ncols-1:
                out_f.write('\t')

            out_f.write('\n')
        
        # block_str = gvc.binarization.simd_tensor_to_txt(allele_tensor, phasing_tensor)
        # out_f.write(block_str)

        log.info(time.time() - stime)


def _get_tensor_shape(enc_var, param_set):

    shapes = []
    for i_bin_mat in range(param_set.num_variants_flags):

        if param_set.variants_coder_ids[i_bin_mat] == 0:
            header_bytes = enc_var.variants_payloads[i_bin_mat].read(jbig.BIE_HEADER_LEN)
            bin_mat_nrows, bin_mat_ncols = jbig.get_shape(header_bytes)

        elif param_set.variants_coder_ids[i_bin_mat] == 1:
            raise NotImplementedError()
            header_bytes = enc_var.variants_payloads[i_bin_mat].read(ds.GabacMatrix.NROW_LEN + ds.GabacMatrix.NCOL_LEN)

            bin_mat_nrows = gvc.utils.bytes_to_int(header_bytes[:4])
            bin_mat_ncols = gvc.utils.bytes_to_int(header_bytes[4:])
        else:
            raise NotImplementedError()

        if param_set.transpose_variants_mat_flags[i_bin_mat]:

            bin_mat_nrows, bin_mat_ncols = bin_mat_ncols, bin_mat_nrows

        shapes.append((bin_mat_nrows, bin_mat_ncols))

    nrows, ncols = shapes[0]
    for shape in shapes[1:]:
        assert nrows == shape[0]
        assert nrows == shape[1]

    if param_set.binarization_flag == 0:

        if param_set.concat_axis == 0:
            nrows = nrows // param_set.num_bin_mat
        elif param_set.concat_axis == 1:
            ncols = ncols // param_set.num_bin_mat

    elif param_set.binarization_flag == 1:
        vector_amax = ds.VectorAMax.from_bytes(enc_var.variants_amax_payload.read()).vector
        nrows += int(np.sum(vector_amax[vector_amax != 1]))

    tensor_nrows = nrows
    tensor_ncols = ncols // param_set.p
    tensor_nchannels = param_set.p

    # assert tensor_ncols * tensor_nchannels == ncols

    return (tensor_nrows, tensor_ncols, tensor_nchannels)

class DecoderContext(object):
    def __init__(self):
        self.parameter_sets = {}
        self.access_units = {}

        self.nrows = []
        self.ncols = None

        self.curr_access_unit = None
        self.curr_parameter_set = None

    def set_access_unit(self, access_unit_id):
        try:
            self.curr_access_unit = self.access_units[access_unit_id]
        except KeyError:
            log.error('No access unit with id {} found'.format(access_unit_id))
            raise gvc.errors.GvcError()

        parameter_set_id = self.curr_access_unit.header.parameter_set_id
        try:
            self.curr_parameter_set = self.parameter_sets[parameter_set_id]
        except KeyError:
            log.error('No parameter set with id {} found'.format(parameter_set_id))
            raise gvc.errors.GvcError()

    def __len__(self):
        return len(self.parameter_sets)

    def __getitem__(self,
        parameter_set_id: int
    ) -> ds.ParameterSet:

        return self.parameter_sets[parameter_set_id]

    def __setitem__(self,
                    parameter_set_id: int,
                    value: dict
                    ):
        self.parameter_sets[parameter_set_id] = value

class Decoder(object):
    def __init__(self,
        input_fpath,
        output_fpath=None,
    ):

        self.input_fpath = input_fpath
        self.output_fpath = output_fpath
        self._f = open(self.input_fpath, 'rb')
        self._bitstream_reader = gvc.bitstream.BitstreamReader(self._f)

        self.decoder_context = DecoderContext()
        self._cache_data()

        self.index = ds.Index.from_gvc_fpath(input_fpath, self.decoder_context)

    def _decode_parameter_set(self):
        log.info('decoding parameter set')

        header = ds.DataUnitHeader.from_bitstream(consts.DataUnitType.PARAMETER_SET,
                                               self._bitstream_reader)
        param_set = ds.ParameterSet.from_bitstream(self._bitstream_reader,
                                                                    header)
        self.decoder_context.parameter_sets[param_set.parameter_set_id] = param_set

    def _cache_access_unit(self):
        start_pos = self._f.tell()

        # data_unit_size = self._bitstream_reader.read_bits(consts.DATA_UNIT_SIZE_LEN * 8)  # data_unit_size
        # data_unit_size = self._bitstream_reader.read_bytes(consts.DATA_UNIT_SIZE_LEN)
        acc_unit = ds.AccessUnit.from_bitstream(
            self._bitstream_reader, self.decoder_context.parameter_sets)

        end_pos = self._f.tell()

        # assert end_pos-start_pos+consts.DATA_UNIT_TYPE_LEN == data_unit_size

        param_set = self.decoder_context[acc_unit.header.parameter_set_id]

        nrows = 0
        for block in acc_unit.blocks:
            enc_var = block.block_payload
            tensor_shape = _get_tensor_shape(enc_var, param_set)

            if tensor_shape[2] != param_set.p:
                raise ValueError()

            nrows += tensor_shape[0]

            if self.decoder_context.ncols is not None:
                if tensor_shape[1] != self.decoder_context.ncols or tensor_shape[2] != param_set.p: 
                    log.error('Tensor shape not consistent')
            else:
                self.decoder_context.ncols = tensor_shape[1]

        #? Preallocate list with values None
        #? Created to handle case where access unit id stored is not in order
        access_unit_id = acc_unit.header.access_unit_id
        if len(self.decoder_context.nrows) <= access_unit_id:
            for _ in range(access_unit_id - len(self.decoder_context.nrows)):
                self.decoder_context.nrows.append(None)

            self.decoder_context.nrows.append(nrows)
        else:
            self.decoder_context.nrows[access_unit_id] = nrows
    
        self.decoder_context.access_units[access_unit_id] = acc_unit
        
        self._f.seek(end_pos)

    def _cache_data(self):
        log.info('Caching data')

        while True:
            data_unit_type = self._bitstream_reader.read_bits(consts.DATA_UNIT_TYPE_LEN * 8)

            if not self._bitstream_reader.read:
                break

            if data_unit_type == consts.DataUnitType.PARAMETER_SET:
                self._decode_parameter_set()

            elif data_unit_type == consts.DataUnitType.ACCESS_UNIT:
                self._cache_access_unit()
                
            else:
                raise TypeError('invalid data unit type: {}'.format(data_unit_type))

            self._bitstream_reader._reset()

        for i_access_unit in range(self.num_access_units):
            try:
                self.decoder_context.access_units[i_access_unit]
            except KeyError:
                log.error('Missing access_unit_id: {}'.format(i_access_unit))
                raise ValueError('Missing access_unit_id: {}'.format(i_access_unit))

    @property
    def num_access_units(self):
        return len(self.decoder_context.access_units)

    @property
    def num_parameter_sets(self):
        return len(self.decoder_context.parameter_sets)

    @property
    def num_decoded_data_units(self):
        return self.num_access_units + self.num_parameter_sets

    def decode_all(self):

        with open(self.output_fpath, 'w') as out_f:

            for i_access_unit in range(self.num_access_units):
                self.decoder_context.set_access_unit(i_access_unit)

                _decode_and_write_access_unit(out_f, self.decoder_context)
                
    def random_access(self, start_pos, end_pos, sample_ids=None):
        
        if end_pos == -1:
            end_pos = start_pos
            
        query_col_ids = self.index.query_columns(sample_ids)

        #? Query block and parameter set id given position
        blk_ps_id_pairs = self.index.query_blk(start_pos, end_pos)

        if blk_ps_id_pairs.shape[0]:

            for i_block in range(blk_ps_id_pairs.shape[0]):
                block_id, block, param_set_id = blk_ps_id_pairs[i_block, :]

                row_slice = self.index.get_row_mask(block_id, start_pos, end_pos)
                if row_slice.start < row_slice.stop:
                    param_set = self.decoder_context.parameter_sets[param_set_id]

                    _selective_decode_encoded_variants(
                        param_set, 
                        block.block_payload, 
                        row_slice, 
                        query_col_ids
                    )

                #? No variant found given POSs
                else:
                    # return None
                    pass

                pass

        #? No block found given POSs
        else:
            return None

    def compare(self):

        with open(self.output_fpath, 'r') as out_f:

            for i_access_unit in range(self.num_access_units):
                self.decoder_context.set_access_unit(i_access_unit)

                allele_tensor, phasing_tensor = _decode_access_unit(self.decoder_context)
                log.info("Shape: {} {}".format(allele_tensor.shape, phasing_tensor.shape))

                gvc.binarization.compare_tensor_to_txt(allele_tensor, phasing_tensor, out_f)

    # def stat(self):
    #     try:
    #         import pandas as pd
    #     except:
    #         raise RuntimeError("Please install pandas!")
        
    #     df = pd.DataFrame(columns=[
    #         "AccessUnitID", "BlockID", 
    #         "NumAlleleBinMat", "AlleleBinMat", "AlleleRowIds", "AlleleColIds", 
    #         "AMax", 
    #         "PhaseBinMat", "PhaseRowIds", "PhaseColIds",
    #     ])

    #     for access_unit_id in range(self.num_access_units):
    #         self.decoder_context.set_access_unit(access_unit_id)

    #         for block_id, block in enumerate(self.decoder_context.curr_access_unit.blocks):
    #             log.info('Stat block {}'.format(block_id))

    #             stat = block.stat()

    #             stat["AccessUnitID"] = access_unit_id
    #             stat["BlockID"] = block_id
                
    #             df = df.append(stat, ignore_index=True)

    #     df.to_csv(
    #         self.output_fpath, index=False
    #     )

    # def cprofile(self):

    #     for i_access_unit in range(self.num_access_units):
    #         self.decoder_context.set_access_unit(i_access_unit)

    #         for i, block in enumerate(self.decoder_context.curr_access_unit.blocks):
    #             log.info('Decoding block {}'.format(i))
    #             allele_tensor, phasing_tensor = _decode_block_payload(self.decoder_context.curr_parameter_set, block)






