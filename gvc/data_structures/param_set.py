from dataclasses import dataclass
import logging as log
import typing as t
from . import consts
from .data_unit import DataUnitHeader
from .. import bitstream

class ParameterSet():
    def __init__(
        self,
        parameter_set_id: int,
        any_missing_flag:bool,
        not_available_flag:bool,
        p:int,
        binarization_id:int, 
        num_bin_mat: int,
        concat_axis: int,
        sort_variants_row_flags: t.List[bool],
        sort_variants_col_flags: t.List[bool],
        transpose_variants_mat_flags: t.List[bool],
        variants_coder_ids: t.List[str],
        encode_phase_data: bool,
        phase_value:bool=None,
        sort_phases_row_flag:bool=None,
        sort_phases_col_flag:bool=None,
        transpose_phase_mat_flag:bool=None,
        phase_coder_ids:str=None
    ):
        self.parameter_set_id = parameter_set_id

        self.any_missing_flag = any_missing_flag
        self.not_available_flag = not_available_flag
        self.p = p
        self.binarization_id = binarization_id

        self.num_bin_mat = None
        self.concat_axis = None
        # Binarization using bit plane
        if self.binarization_id in [consts.BinarizationID.BIT_PLANE]:
            self.num_bin_mat = num_bin_mat
            # 2 bits : 0, 1, 2 (do not concatenate)
            if concat_axis not in (0, 1, 2):
                log.error('Invalid value for concat_axis')
                raise ValueError('Invalid value for concat_axis')

            self.concat_axis = concat_axis

            if self.concat_axis in (0, 1):
                self.num_variants_flags = 1

            elif self.concat_axis == 2:
                self.num_variants_flags = num_bin_mat         

        # Binarization by row splitting
        elif self.binarization_id in [consts.BinarizationID.ROW_BIN_SPLIT]:
            self.num_variants_flags = 1
        else:
            raise ValueError('Invalid binarization flag')

        if len(sort_variants_row_flags) == self.num_variants_flags and len(sort_variants_col_flags) == self.num_variants_flags and \
            len(transpose_variants_mat_flags) == self.num_variants_flags and len(variants_coder_ids) == self.num_variants_flags:

            self.sort_variants_row_flags = sort_variants_row_flags
            self.sort_variants_col_flags = sort_variants_col_flags
            self.transpose_variants_mat_flags = transpose_variants_mat_flags
            self.variants_coder_ids = variants_coder_ids

        else:
            log.error('Invalid variants flags')
            raise ValueError('Invalid variants flags')

        self.encode_phase_data = encode_phase_data

        # Initialize default value of phase matrix flags
        self.phase_value = None
        self.sort_phases_row_flag = None
        self.sort_phases_col_flag = None
        self.transpose_phase_mat_flag = None

        if self.encode_phase_data:
            self.sort_phases_row_flag = sort_phases_row_flag
            self.sort_phases_col_flag = sort_phases_col_flag
            self.transpose_phase_mat_flag = transpose_phase_mat_flag
            self.phase_coder_ids = phase_coder_ids
        else:
            self.phase_value = phase_value

    def __eq__(self, 
        other
    ):
        """

        Parameters
        ----------
        other: gvc.data_structures.ParameterSet

        Returns
        -------

        """

        if other is None:
            return False

        if self.any_missing_flag != other.any_missing_flag:
            return False
        if self.not_available_flag != other.not_available_flag:
            return False
        if self.p != other.p:
            return False
        if self.binarization_id != other.binarization_id:
            return False
        if self.num_bin_mat != other.num_bin_mat:
            return False
        if self.concat_axis != other.concat_axis:
            return False
        if self.num_variants_flags != other.num_variants_flags:
            return False

        for i in range(self.num_variants_flags):
            if self.sort_variants_row_flags[i] != other.sort_variants_row_flags[i]:
                return False
            if self.sort_variants_col_flags[i] != other.sort_variants_col_flags[i]:
                return False
            if self.transpose_variants_mat_flags[i] != other.transpose_variants_mat_flags[i]:
                return False
            if self.variants_coder_ids[i] != other.variants_coder_ids[i]:
                return False

        if self.encode_phase_data == other.encode_phase_data:
            if self.encode_phase_data:
                if self.sort_phases_col_flag != other.sort_phases_col_flag:
                    return False
                if self.sort_phases_row_flag != other.sort_phases_row_flag:
                    return False
                if self.transpose_phase_mat_flag != other.transpose_phase_mat_flag:
                    return False
                if self.phase_coder_ids != other.phase_coder_ids:
                    return False
            else:
                if self.phase_value != other.phase_value:
                    return False
        else:
            return False            

        return True

    def to_bitio(self):
        """

        Parameters
        ----------
        with_header: bool, optional, default: True
            Generate header binary such as DATA_UNIT_TYPE and length of data unit

        Returns
        -------
        header_bitio : gvc.bitstream.BitIO
            BitIO object of header (DATA_UNIT_TYPE and total length of data unit)
        data_bitio : gvc.bitstream.BitIO
            BitIO object of the data of parameter set
        """

        data_bitio = bitstream.BitIO()

        data_bitio.write(self.parameter_set_id, consts.PARAMETER_SET_ID_LEN*8)

        data_bitio.write(self.any_missing_flag, consts.ANY_MISSING_FLAG_BITLEN)
        data_bitio.write(self.not_available_flag, consts.NOT_AVAILABLE_FLAG_BITLEN)

        # Minimum value of p is 1, thus shifted by 1 to reduce bits required
        data_bitio.write(self.p-1, consts.P_BITLEN)

        data_bitio.write(self.binarization_id, consts.BINARIZAION_ID_BITLEN)

        if self.binarization_id == 0:
            data_bitio.write(self.num_bin_mat, consts.NUM_BIN_MAT_BITLEN)
            data_bitio.write(self.concat_axis, consts.CONCAT_AXIS_BITLEN)

        data_bitio.write(self.encode_phase_data, consts.ENCODE_PHASE_DATA_BITLEN)

        for i in range(self.num_variants_flags):
            data_bitio.write(self.sort_variants_row_flags[i], consts.SORT_VARIANTS_FLAG_BITLEN)
            data_bitio.write(self.sort_variants_col_flags[i], consts.SORT_VARIANTS_FLAG_BITLEN)
            data_bitio.write(self.transpose_variants_mat_flags[i], consts.TRANSPOSE_FLAG_BITLEN)
            data_bitio.write(self.variants_coder_ids[i], consts.CODER_ID_BITLEN)

        if self.encode_phase_data:
            data_bitio.write(self.sort_phases_row_flag, consts.SORT_VARIANTS_FLAG_BITLEN)
            data_bitio.write(self.sort_phases_col_flag, consts.SORT_VARIANTS_FLAG_BITLEN)
            data_bitio.write(self.transpose_phase_mat_flag, consts.TRANSPOSE_FLAG_BITLEN)
            data_bitio.write(self.phase_coder_ids, consts.CODER_ID_BITLEN)
        else:
            data_bitio.write(self.phase_value, consts.PHASE_VALUE_BITLEN)

        data_bitio.align_to_byte()

        return data_bitio

    def to_bytes(self, header=True):
        # header_bitio, data_bitio = self.to_bitio()
        payload_b = self.to_bitio().to_bytes()

        if header:
            header = DataUnitHeader(consts.DataUnitType.PARAMETER_SET, len(payload_b))

            return bytes(header.to_barray()) + payload_b
        else:
            return payload_b

    def to_bitstream(self, bitstream_writer):
        raise NotImplementedError()

    @classmethod
    def from_bitstream(
        cls, 
        bitstream_reader:bitstream.BitstreamReader,
        header:DataUnitHeader,
    ):
        """

        Parameters
        ----------

        Returns
        -------

        """

        start_pos = bitstream_reader.tell()

        parameter_set_id = bitstream_reader.read_bytes(consts.PARAMETER_SET_ID_LEN, ret_int=True)

        any_missing_flag = bitstream_reader.read_bits(consts.ANY_MISSING_FLAG_BITLEN)
        not_available_flag = bitstream_reader.read_bits(consts.NOT_AVAILABLE_FLAG_BITLEN)

        # Revert shifted value
        p = bitstream_reader.read_bits(consts.P_BITLEN) + 1

        binarization_id = bitstream_reader.read_bits(consts.BINARIZAION_ID_BITLEN)

        num_bin_mat = None
        concat_axis = None
        if binarization_id == consts.BinarizationID.BIT_PLANE:
            num_bin_mat = bitstream_reader.read_bits(consts.NUM_BIN_MAT_BITLEN)
            concat_axis = bitstream_reader.read_bits(consts.CONCAT_AXIS_BITLEN)            

        encode_phase_data = bitstream_reader.read_bits(consts.ENCODE_PHASE_DATA_BITLEN)

        # Binarization using bit plane
        if binarization_id == consts.BinarizationID.BIT_PLANE:
            if concat_axis in (0, 1):
                num_variants_flags = 1

            elif concat_axis == 2:
                num_variants_flags = num_bin_mat         

        # Binarization by row splitting
        elif binarization_id in [consts.BinarizationID.ROW_BIN_SPLIT]:
            num_variants_flags = 1
        
        sort_variants_row_flags = []
        sort_variants_col_flags = []
        transpose_variants_mat_flags = []
        variants_coder_ids = []

        for _ in range(num_variants_flags):
            sort_variants_row_flags.append(bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN))
            sort_variants_col_flags.append(bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN))
            transpose_variants_mat_flags.append(bitstream_reader.read_bits(consts.TRANSPOSE_FLAG_BITLEN))
            variants_coder_ids.append(bitstream_reader.read_bits(consts.CODER_ID_BITLEN))

        phase_value = None
        sort_phases_row_flag = None
        sort_phases_col_flag = None
        transpose_phase_mat_flag = None
        phase_coder_ids = None

        if encode_phase_data:
            sort_phases_row_flag = bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN)
            sort_phases_col_flag = bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN)
            transpose_phase_mat_flag = bitstream_reader.read_bits(consts.TRANSPOSE_FLAG_BITLEN)
            phase_coder_ids = bitstream_reader.read_bits(consts.CODER_ID_BITLEN)
        else:
            phase_value = bitstream_reader.read_bits(consts.PHASE_VALUE_BITLEN)

        bitstream_reader.align_to_byte()

        #? DATA_UNIT read outside this function, thus offset DATA_UNIT_TYPE_LEN is required
        # assert (bitstream_reader.tell() - start_pos + consts.DATA_UNIT_TYPE_LEN) == data_len
        assert (bitstream_reader.tell() - start_pos) == header.content_len

        return cls(
            parameter_set_id,
            any_missing_flag,
            not_available_flag,
            p,
            binarization_id,
            num_bin_mat,
            concat_axis,
            sort_variants_row_flags,
            sort_variants_col_flags,
            transpose_variants_mat_flags,
            variants_coder_ids,
            encode_phase_data,
            phase_value,
            sort_phases_row_flag,
            sort_phases_col_flag,
            transpose_phase_mat_flag,
            phase_coder_ids
        )