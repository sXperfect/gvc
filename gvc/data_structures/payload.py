import logging as log
import typing as t
from . import consts
from .param_set import ParameterSet
from ..bitstream import RandomAccessHandler
from .. import utils

class GenotypePayload(object):
    def __init__(self,
        param_set:ParameterSet,
        variants_payloads:t.List[bytes],
        variants_row_ids_payloads:t.List[bytes],
        variants_col_ids_payloads:t.List[bytes],
        variants_amax_payload=None,
        phase_payload:bytes=None,
        phase_row_ids_payload:bytes=None, 
        phase_col_ids_payload:bytes=None,
        missing_rep_val:int=None,
        na_rep_val:int=None,
    ):

        if not isinstance(param_set, ParameterSet):
            raise TypeError('param_set is not a ParameterSet object')

        if param_set.num_variants_flags != len(variants_payloads) or param_set.num_variants_flags != len(variants_row_ids_payloads) \
            or param_set.num_variants_flags != len(variants_col_ids_payloads):
            raise ValueError('Inconsistency detected between flags and payloads')

        self.variants_payloads = variants_payloads
        self.variants_row_ids_payloads = variants_row_ids_payloads
        self.variants_col_ids_payloads = variants_col_ids_payloads
        if param_set.binarization_id in (0,3):
            self.variants_amax_payload = None
        elif param_set.binarization_id in (1,2) and variants_amax_payload is not None:
            self.variants_amax_payload = variants_amax_payload
        else:
            raise ValueError('Invalid binarization flag combination')

        self.phase_payload = None
        self.phase_row_ids_payload = None
        self.phase_col_ids_payload = None

        if phase_payload is not None:
            self.phase_payload = phase_payload
            self.phase_row_ids_payload = phase_row_ids_payload
            self.phase_col_ids_payload = phase_col_ids_payload
            

        self.missing_rep_val = missing_rep_val
        self.na_rep_val = na_rep_val

    def __len__(self):

        length = 0

        num_variants_flags = len(self.variants_payloads)
        for i_bin_mat in range(num_variants_flags):
            
            length += consts.VARIANTS_PAYLOAD_SIZES_LEN
            length += len(self.variants_payloads[i_bin_mat])

            # Concate bytes if value is not None
            if self.variants_row_ids_payloads[i_bin_mat] is not None:
                length += consts.ROW_IDS_SIZE_LEN
                length += len(self.variants_row_ids_payloads[i_bin_mat])

            # Concate bytes if value is not None
            if self.variants_col_ids_payloads[i_bin_mat] is not None:
                length += consts.COL_IDS_SIZE_LEN
                length += len(self.variants_col_ids_payloads[i_bin_mat])

        if self.variants_amax_payload is not None:
            length += consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN
            length += len(self.variants_amax_payload)

        if self.phase_payload is not None:

            length += consts.PHASE_PAYLOAD_SIZE_LEN
            length += len(self.phase_payload)

            # Concate bytes if value is not None
            if self.phase_row_ids_payload:
                length += consts.ROW_IDS_SIZE_LEN
                length += len(self.phase_row_ids_payload)

            # Concate bytes if value is not None
            if self.phase_col_ids_payload:
                length += consts.COL_IDS_SIZE_LEN
                length += len(self.phase_col_ids_payload)
                
        if self.missing_rep_val is not None:
            length += 1
            
        if self.na_rep_val is not None:
            length += 1

        return length

    def stat(self):

        total_allele_payload_size = 0
        total_allele_row_ids_payload_size = 0
        total_allele_col_ids_payload_size = 0

        num_variants_flags = len(self.variants_payloads)
        for i_bin_mat in range(num_variants_flags):
            
            total_allele_payload_size += consts.VARIANTS_PAYLOAD_SIZES_LEN
            total_allele_payload_size += len(self.variants_payloads[i_bin_mat])

            # Concate bytes if value is not None
            if self.variants_row_ids_payloads[i_bin_mat] is not None:
                total_allele_row_ids_payload_size += consts.ROW_IDS_SIZE_LEN
                total_allele_row_ids_payload_size += len(self.variants_row_ids_payloads[i_bin_mat])

            # Concate bytes if value is not None
            if self.variants_col_ids_payloads[i_bin_mat] is not None:
                total_allele_col_ids_payload_size += consts.COL_IDS_SIZE_LEN
                total_allele_col_ids_payload_size += len(self.variants_col_ids_payloads[i_bin_mat])

        total_amax_payload_size = 0
        if self.variants_amax_payload is not None:
            total_amax_payload_size += consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN
            total_amax_payload_size += len(self.variants_amax_payload)

        total_phase_payload_size = 0
        total_phase_row_ids_payload_size = 0
        total_phase_col_ids_payload_size = 0

        if self.phase_payload is not None:

            total_phase_payload_size += consts.PHASE_PAYLOAD_SIZE_LEN
            total_phase_payload_size += len(self.phase_payload)

            # Concate bytes if value is not None
            if self.phase_row_ids_payload:
                total_phase_row_ids_payload_size += consts.ROW_IDS_SIZE_LEN
                total_phase_row_ids_payload_size += len(self.phase_row_ids_payload)

            # Concate bytes if value is not None
            if self.phase_col_ids_payload:
                total_allele_col_ids_payload_size += consts.COL_IDS_SIZE_LEN
                total_allele_col_ids_payload_size += len(self.phase_col_ids_payload)

        #TODO: Both missing_rep_val and na_rep_val are not yet included
        return {
            "NumAlleleBinMat" : num_variants_flags,
            "AlleleBinMat": total_allele_payload_size,
            "AlleleRowIds": total_allele_row_ids_payload_size,
            "AlleleColIds": total_allele_col_ids_payload_size,
            "AMax": total_amax_payload_size,
            "PhaseBinMat": total_phase_payload_size,
            "PhaseRowIds": total_phase_row_ids_payload_size,
            "PhaseColIds": total_phase_col_ids_payload_size,
        }

    def to_barray(self):
        payload = bytearray()

        num_variants_flags = len(self.variants_payloads)
        for i_bin_mat in range(num_variants_flags):
            
            payload += utils.int2bstr(len(self.variants_payloads[i_bin_mat]), consts.VARIANTS_PAYLOAD_SIZES_LEN)
            payload += self.variants_payloads[i_bin_mat]

            # Concate bytes if value is not None
            if self.variants_row_ids_payloads[i_bin_mat] is not None:
                payload += utils.int2bstr(len(self.variants_row_ids_payloads[i_bin_mat]), consts.ROW_IDS_SIZE_LEN)
                payload += self.variants_row_ids_payloads[i_bin_mat]

            # Concate bytes if value is not None
            if self.variants_col_ids_payloads[i_bin_mat]:
                payload += utils.int2bstr(len(self.variants_col_ids_payloads[i_bin_mat]), consts.COL_IDS_SIZE_LEN)
                payload += self.variants_col_ids_payloads[i_bin_mat]

        if self.variants_amax_payload is not None:
            payload += utils.int2bstr(len(self.variants_amax_payload), consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN)
            payload += self.variants_amax_payload

        if self.phase_payload is not None:

            payload += utils.int2bstr(len(self.phase_payload), consts.PHASE_PAYLOAD_SIZE_LEN)
            payload += self.phase_payload

            # Concate bytes if value is not None
            if self.phase_row_ids_payload is not None:
                payload += utils.int2bstr(len(self.phase_row_ids_payload), consts.ROW_IDS_SIZE_LEN)
                payload += self.phase_row_ids_payload

            # Concate bytes if value is not None
            if self.phase_col_ids_payload is not None:
                payload += utils.int2bstr(len(self.phase_col_ids_payload), consts.COL_IDS_SIZE_LEN)
                payload += self.phase_col_ids_payload
                
        if self.missing_rep_val is not None:
            payload += utils.int2bstr(self.missing_rep_val, consts.MISSING_REP_VAL_LEN)

        if self.na_rep_val is not None:
            payload += utils.int2bstr(self.na_rep_val, consts.NA_REP_VAL_LEN)

        return payload

    def to_bytes(self):
        return bytes(self.to_barray())

    @classmethod
    def from_bitstream(
        cls, 
        reader,
        param_set:ParameterSet, 
        block_payload_size
    ):

        start_pos = reader.tell()

        variants_payloads = []
        variants_row_ids_payloads = []
        variants_col_ids_payloads = []

        for i_bin_mat in range(param_set.num_variants_flags):
            variants_payload_size = reader.read_bytes(consts.VARIANTS_PAYLOAD_SIZES_LEN, ret_int=True)

            variants_payload = RandomAccessHandler(
                reader,
                reader.tell(),
                variants_payload_size
            )
            reader.seek(variants_payload_size, 1)

            variants_payloads.append(variants_payload)

            variants_row_ids_payload = None
            if param_set.sort_variants_row_flags[i_bin_mat]:
                
                variants_row_ids_payload_size = reader.read_bytes(consts.ROW_IDS_SIZE_LEN, ret_int=True)

                variants_row_ids_payload = RandomAccessHandler(
                    reader,
                    reader.tell(),
                    variants_row_ids_payload_size
                )
                reader.seek(variants_row_ids_payload_size, 1)

            variants_row_ids_payloads.append(variants_row_ids_payload)

            variants_col_ids_payload = None
            if param_set.sort_variants_col_flags[i_bin_mat]:
                variants_col_ids_payload_size = reader.read_bytes(consts.COL_IDS_SIZE_LEN, ret_int=True)
                

                variants_col_ids_payload = RandomAccessHandler(
                    reader,
                    reader.tell(),
                    variants_col_ids_payload_size
                )
                reader.seek(variants_col_ids_payload_size, 1)

            variants_col_ids_payloads.append(variants_col_ids_payload)

        variants_amax_payload = None
        if param_set.binarization_id in [1, 2]:
            variants_amax_payload_size = reader.read_bytes(consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN, ret_int=True)

            variants_amax_payload = RandomAccessHandler(
                reader,
                reader.tell(),
                variants_amax_payload_size
            )
            reader.seek(variants_amax_payload_size, 1)

        phase_payload = None
        phase_row_ids_payload = None
        phase_col_ids_payload = None

        if param_set.encode_phase_data:
            phase_payload_size = reader.read_bytes(consts.PHASE_PAYLOAD_SIZE_LEN, ret_int=True)

            phase_payload = RandomAccessHandler(
                reader,
                reader.tell(),
                phase_payload_size
            )
            reader.seek(phase_payload_size, 1)

            if param_set.sort_phases_row_flag:
                phase_row_ids_payload_size = reader.read_bytes(consts.ROW_IDS_SIZE_LEN, ret_int=True)

                phase_row_ids_payload = RandomAccessHandler(
                    reader,
                    reader.tell(),
                    phase_row_ids_payload_size
                )
                reader.seek(phase_row_ids_payload_size, 1)

            if param_set.sort_phases_col_flag:
                phase_col_ids_payload_size = reader.read_bytes(consts.COL_IDS_SIZE_LEN, ret_int=True)
                
                phase_col_ids_payload = RandomAccessHandler(
                    reader,
                    reader.tell(),
                    phase_col_ids_payload_size
                )
                reader.seek(phase_col_ids_payload_size, 1)
                
        if param_set.any_missing_flag:
            missing_rep_val = reader.read_bytes(consts.MISSING_REP_VAL_LEN, ret_int=True)
        else:
            missing_rep_val = None
            
        if param_set.not_available_flag:
            na_rep_val = reader.read_bytes(consts.NA_REP_VAL_LEN, ret_int=True)
        else:
            na_rep_val = None

        end_pos = reader.tell()

        enc_var =  cls(
            param_set,
            variants_payloads,
            variants_row_ids_payloads=variants_row_ids_payloads,
            variants_col_ids_payloads=variants_col_ids_payloads,
            variants_amax_payload=variants_amax_payload,
            phase_payload=phase_payload,
            phase_row_ids_payload=phase_row_ids_payload,
            phase_col_ids_payload=phase_col_ids_payload,
            missing_rep_val=missing_rep_val,
            na_rep_val=na_rep_val,
        )

        assert len(enc_var) == block_payload_size and (end_pos-start_pos) == block_payload_size
        
        return enc_var