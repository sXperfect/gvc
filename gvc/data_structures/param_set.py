from . import consts
from .data_unit import DataUnitHeader
from .. import bitstream


def _fits_unsigned(value, bits):
    return isinstance(value, int) and not isinstance(value, bool) and 0 <= value < (1 << bits)


def _flag(value, name):
    if value not in (0, 1, False, True):
        raise ValueError("{} must be 0 or 1".format(name))
    return bool(value)


class ParameterSet:
    def __init__(
        self,
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
        phase_value=None,
        sort_phases_row_flag=None,
        sort_phases_col_flag=None,
        transpose_phase_mat_flag=None,
        phase_coder_ids=None,
    ):
        if not _fits_unsigned(parameter_set_id, consts.PARAMETER_SET_ID_LEN * 8):
            raise ValueError("parameter_set_id is outside the serialized range")
        if not isinstance(p, int) or isinstance(p, bool) or not 1 <= p <= (1 << consts.P_BITLEN):
            raise ValueError("p must be in the serialized range 1..{}".format(1 << consts.P_BITLEN))

        try:
            binarization_id = consts.BinarizationID(binarization_id)
        except ValueError as exc:
            raise ValueError("invalid binarization id: {}".format(binarization_id)) from exc

        self.parameter_set_id = parameter_set_id
        self.any_missing_flag = _flag(any_missing_flag, "any_missing_flag")
        self.not_available_flag = _flag(not_available_flag, "not_available_flag")
        self.p = p
        self.binarization_id = binarization_id

        self.num_bin_mat = None
        self.concat_axis = None

        if binarization_id == consts.BinarizationID.BIT_PLANE:
            if (
                not isinstance(num_bin_mat, int)
                or isinstance(num_bin_mat, bool)
                or not 1 <= num_bin_mat < (1 << consts.NUM_BIN_MAT_BITLEN)
            ):
                raise ValueError("num_bin_mat is outside the serialized range")
            if concat_axis not in (0, 1, 2):
                raise ValueError("concat_axis must be 0, 1, or 2")

            self.num_bin_mat = num_bin_mat
            self.concat_axis = int(concat_axis)
            self.num_variants_flags = 1 if concat_axis in (0, 1) else num_bin_mat

        elif binarization_id == consts.BinarizationID.ROW_BIN_SPLIT:
            self.num_variants_flags = 1

        flag_lengths = {
            len(sort_variants_row_flags),
            len(sort_variants_col_flags),
            len(transpose_variants_mat_flags),
            len(variants_coder_ids),
        }
        if flag_lengths != {self.num_variants_flags}:
            raise ValueError("variant flag arrays do not match num_variants_flags")

        coder_limit = 1 << consts.CODER_ID_BITLEN
        if any(
            not isinstance(coder_id, int)
            or isinstance(coder_id, bool)
            or not 0 <= coder_id < coder_limit
            for coder_id in variants_coder_ids
        ):
            raise ValueError("variant coder id is outside the serialized range")

        self.sort_variants_row_flags = [
            _flag(v, "sort_variants_row_flags") for v in sort_variants_row_flags
        ]
        self.sort_variants_col_flags = [
            _flag(v, "sort_variants_col_flags") for v in sort_variants_col_flags
        ]
        self.transpose_variants_mat_flags = [
            _flag(v, "transpose_variants_mat_flags")
            for v in transpose_variants_mat_flags
        ]
        self.variants_coder_ids = list(variants_coder_ids)

        self.encode_phase_data = _flag(encode_phase_data, "encode_phase_data")
        self.phase_value = None
        self.sort_phases_row_flag = None
        self.sort_phases_col_flag = None
        self.transpose_phase_mat_flag = None
        self.phase_coder_ids = None

        if self.encode_phase_data:
            if None in (
                sort_phases_row_flag,
                sort_phases_col_flag,
                transpose_phase_mat_flag,
                phase_coder_ids,
            ):
                raise ValueError("phase coding flags must be provided when phase data is encoded")
            if (
                not isinstance(phase_coder_ids, int)
                or isinstance(phase_coder_ids, bool)
                or not 0 <= phase_coder_ids < coder_limit
            ):
                raise ValueError("phase coder id is outside the serialized range")

            self.sort_phases_row_flag = _flag(
                sort_phases_row_flag, "sort_phases_row_flag"
            )
            self.sort_phases_col_flag = _flag(
                sort_phases_col_flag, "sort_phases_col_flag"
            )
            self.transpose_phase_mat_flag = _flag(
                transpose_phase_mat_flag, "transpose_phase_mat_flag"
            )
            self.phase_coder_ids = phase_coder_ids
        else:
            if phase_value not in (0, 1, False, True):
                raise ValueError("phase_value must be 0 or 1 when phase data is not encoded")
            self.phase_value = _flag(phase_value, "phase_value")

    def __eq__(self, other):
        if not isinstance(other, ParameterSet):
            return False

        semantic_fields = (
            "any_missing_flag",
            "not_available_flag",
            "p",
            "binarization_id",
            "num_bin_mat",
            "concat_axis",
            "num_variants_flags",
            "sort_variants_row_flags",
            "sort_variants_col_flags",
            "transpose_variants_mat_flags",
            "variants_coder_ids",
            "encode_phase_data",
            "phase_value",
            "sort_phases_row_flag",
            "sort_phases_col_flag",
            "transpose_phase_mat_flag",
            "phase_coder_ids",
        )
        return all(getattr(self, name) == getattr(other, name) for name in semantic_fields)

    def to_bitio(self):
        data_bitio = bitstream.BitIO()
        data_bitio.write(self.parameter_set_id, consts.PARAMETER_SET_ID_LEN * 8)
        data_bitio.write(self.any_missing_flag, consts.ANY_MISSING_FLAG_BITLEN)
        data_bitio.write(self.not_available_flag, consts.NOT_AVAILABLE_FLAG_BITLEN)
        data_bitio.write(self.p - 1, consts.P_BITLEN)
        data_bitio.write(int(self.binarization_id), consts.BINARIZAION_ID_BITLEN)

        if self.binarization_id == consts.BinarizationID.BIT_PLANE:
            data_bitio.write(self.num_bin_mat, consts.NUM_BIN_MAT_BITLEN)
            data_bitio.write(self.concat_axis, consts.CONCAT_AXIS_BITLEN)

        data_bitio.write(self.encode_phase_data, consts.ENCODE_PHASE_DATA_BITLEN)

        for i in range(self.num_variants_flags):
            data_bitio.write(
                self.sort_variants_row_flags[i], consts.SORT_VARIANTS_FLAG_BITLEN
            )
            data_bitio.write(
                self.sort_variants_col_flags[i], consts.SORT_VARIANTS_FLAG_BITLEN
            )
            data_bitio.write(
                self.transpose_variants_mat_flags[i], consts.TRANSPOSE_FLAG_BITLEN
            )
            data_bitio.write(self.variants_coder_ids[i], consts.CODER_ID_BITLEN)

        if self.encode_phase_data:
            data_bitio.write(
                self.sort_phases_row_flag, consts.SORT_VARIANTS_FLAG_BITLEN
            )
            data_bitio.write(
                self.sort_phases_col_flag, consts.SORT_VARIANTS_FLAG_BITLEN
            )
            data_bitio.write(
                self.transpose_phase_mat_flag, consts.TRANSPOSE_FLAG_BITLEN
            )
            data_bitio.write(self.phase_coder_ids, consts.CODER_ID_BITLEN)
        else:
            data_bitio.write(self.phase_value, consts.PHASE_VALUE_BITLEN)

        data_bitio.align_to_byte()
        return data_bitio

    def to_bytes(self, header=True):
        payload = self.to_bitio().to_bytes()
        if not header:
            return payload

        header_obj = DataUnitHeader(consts.DataUnitType.PARAMETER_SET, len(payload))
        return bytes(header_obj.to_barray()) + payload

    @classmethod
    def from_bitstream(cls, bitstream_reader, header):
        start_pos = bitstream_reader.tell()

        parameter_set_id = bitstream_reader.read_bytes(
            consts.PARAMETER_SET_ID_LEN, ret_int=True
        )
        any_missing_flag = bitstream_reader.read_bits(consts.ANY_MISSING_FLAG_BITLEN)
        not_available_flag = bitstream_reader.read_bits(
            consts.NOT_AVAILABLE_FLAG_BITLEN
        )
        p = bitstream_reader.read_bits(consts.P_BITLEN) + 1
        binarization_id = bitstream_reader.read_bits(consts.BINARIZAION_ID_BITLEN)

        num_bin_mat = None
        concat_axis = None
        if binarization_id == consts.BinarizationID.BIT_PLANE:
            num_bin_mat = bitstream_reader.read_bits(consts.NUM_BIN_MAT_BITLEN)
            concat_axis = bitstream_reader.read_bits(consts.CONCAT_AXIS_BITLEN)
            if concat_axis not in (0, 1, 2):
                raise ValueError("invalid serialized concat_axis")
            num_variants_flags = 1 if concat_axis in (0, 1) else num_bin_mat
        elif binarization_id == consts.BinarizationID.ROW_BIN_SPLIT:
            num_variants_flags = 1
        else:
            raise ValueError("unsupported serialized binarization id: {}".format(binarization_id))

        encode_phase_data = bitstream_reader.read_bits(
            consts.ENCODE_PHASE_DATA_BITLEN
        )

        sort_variants_row_flags = []
        sort_variants_col_flags = []
        transpose_variants_mat_flags = []
        variants_coder_ids = []

        for _ in range(num_variants_flags):
            sort_variants_row_flags.append(
                bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN)
            )
            sort_variants_col_flags.append(
                bitstream_reader.read_bits(consts.SORT_VARIANTS_FLAG_BITLEN)
            )
            transpose_variants_mat_flags.append(
                bitstream_reader.read_bits(consts.TRANSPOSE_FLAG_BITLEN)
            )
            variants_coder_ids.append(
                bitstream_reader.read_bits(consts.CODER_ID_BITLEN)
            )

        phase_value = None
        sort_phases_row_flag = None
        sort_phases_col_flag = None
        transpose_phase_mat_flag = None
        phase_coder_ids = None

        if encode_phase_data:
            sort_phases_row_flag = bitstream_reader.read_bits(
                consts.SORT_VARIANTS_FLAG_BITLEN
            )
            sort_phases_col_flag = bitstream_reader.read_bits(
                consts.SORT_VARIANTS_FLAG_BITLEN
            )
            transpose_phase_mat_flag = bitstream_reader.read_bits(
                consts.TRANSPOSE_FLAG_BITLEN
            )
            phase_coder_ids = bitstream_reader.read_bits(consts.CODER_ID_BITLEN)
        else:
            phase_value = bitstream_reader.read_bits(consts.PHASE_VALUE_BITLEN)

        bitstream_reader.align_to_byte()
        consumed = bitstream_reader.tell() - start_pos
        if consumed != header.content_len:
            raise ValueError(
                "parameter-set length mismatch: expected {}, consumed {}".format(
                    header.content_len, consumed
                )
            )

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
            phase_coder_ids,
        )
