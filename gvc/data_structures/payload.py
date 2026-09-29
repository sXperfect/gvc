from numbers import Integral

from . import consts
from .param_set import ParameterSet
from ..bitstream import RandomAccessHandler
from .. import utils


def _payload_bytes(payload):
    if hasattr(payload, "read"):
        return payload.read()
    return bytes(payload)


def _validate_payload_size(payload, size_len, label):
    try:
        size = len(payload)
    except TypeError as exc:
        raise TypeError("{} must have a byte length".format(label)) from exc
    if size < 0 or size >= (1 << (8 * size_len)):
        raise ValueError("{} is outside the serialized size range".format(label))
    return size


def _validate_optional_payload_size(payload, size_len, label):
    if payload is not None:
        _validate_payload_size(payload, size_len, label)


class GenotypePayload:
    def __init__(
        self,
        param_set,
        variants_payloads,
        variants_row_ids_payloads,
        variants_col_ids_payloads,
        variants_amax_payload=None,
        phase_payload=None,
        phase_row_ids_payload=None,
        phase_col_ids_payload=None,
        missing_rep_val=None,
        na_rep_val=None,
    ):
        if not isinstance(param_set, ParameterSet):
            raise TypeError("param_set must be a ParameterSet")

        expected = param_set.num_variants_flags
        if not (
            len(variants_payloads)
            == len(variants_row_ids_payloads)
            == len(variants_col_ids_payloads)
            == expected
        ):
            raise ValueError("variant payload arrays do not match num_variants_flags")

        for i in range(expected):
            if variants_payloads[i] is None:
                raise ValueError("variant payload cannot be None")

            row_expected = param_set.sort_variants_row_flags[i]
            col_expected = param_set.sort_variants_col_flags[i]
            if row_expected != (variants_row_ids_payloads[i] is not None):
                raise ValueError("row permutation payload does not match parameter flags")
            if col_expected != (variants_col_ids_payloads[i] is not None):
                raise ValueError("column permutation payload does not match parameter flags")

        self.variants_payloads = list(variants_payloads)
        self.variants_row_ids_payloads = list(variants_row_ids_payloads)
        self.variants_col_ids_payloads = list(variants_col_ids_payloads)

        for i in range(expected):
            _validate_payload_size(
                self.variants_payloads[i],
                consts.VARIANTS_PAYLOAD_SIZES_LEN,
                "variant payload {}".format(i),
            )
            _validate_optional_payload_size(
                self.variants_row_ids_payloads[i],
                consts.ROW_IDS_SIZE_LEN,
                "variant row permutation {}".format(i),
            )
            _validate_optional_payload_size(
                self.variants_col_ids_payloads[i],
                consts.COL_IDS_SIZE_LEN,
                "variant column permutation {}".format(i),
            )

        if param_set.binarization_id == consts.BinarizationID.BIT_PLANE:
            if variants_amax_payload is not None:
                raise ValueError("bit-plane payload must not contain an AMax vector")
            self.variants_amax_payload = None
        elif param_set.binarization_id == consts.BinarizationID.ROW_BIN_SPLIT:
            if variants_amax_payload is None:
                raise ValueError("row-bin-split payload requires an AMax vector")
            _validate_payload_size(
                variants_amax_payload,
                consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN,
                "AMax payload",
            )
            self.variants_amax_payload = variants_amax_payload
        else:
            raise ValueError("unsupported binarization id")

        if param_set.encode_phase_data:
            if phase_payload is None:
                raise ValueError("phase payload is required by the parameter set")
            if param_set.sort_phases_row_flag != (phase_row_ids_payload is not None):
                raise ValueError("phase row permutation does not match parameter flags")
            if param_set.sort_phases_col_flag != (phase_col_ids_payload is not None):
                raise ValueError("phase column permutation does not match parameter flags")

            _validate_payload_size(
                phase_payload, consts.PHASE_PAYLOAD_SIZE_LEN, "phase payload"
            )
            _validate_optional_payload_size(
                phase_row_ids_payload,
                consts.ROW_IDS_SIZE_LEN,
                "phase row permutation",
            )
            _validate_optional_payload_size(
                phase_col_ids_payload,
                consts.COL_IDS_SIZE_LEN,
                "phase column permutation",
            )
            self.phase_payload = phase_payload
            self.phase_row_ids_payload = phase_row_ids_payload
            self.phase_col_ids_payload = phase_col_ids_payload
        else:
            if any(
                item is not None
                for item in (phase_payload, phase_row_ids_payload, phase_col_ids_payload)
            ):
                raise ValueError("phase payload supplied when phase data is not encoded")
            self.phase_payload = None
            self.phase_row_ids_payload = None
            self.phase_col_ids_payload = None

        if param_set.any_missing_flag != (missing_rep_val is not None):
            raise ValueError("missing-value payload does not match parameter flags")
        if param_set.not_available_flag != (na_rep_val is not None):
            raise ValueError("not-available payload does not match parameter flags")

        for name, value, length in (
            ("missing_rep_val", missing_rep_val, consts.MISSING_REP_VAL_LEN),
            ("na_rep_val", na_rep_val, consts.NA_REP_VAL_LEN),
        ):
            if value is not None:
                if not isinstance(value, Integral) or isinstance(value, bool):
                    raise TypeError("{} must be an integer".format(name))
                value = int(value)
                if not 0 <= value < (1 << (8 * length)):
                    raise ValueError("{} is outside the serialized range".format(name))
                if name == "missing_rep_val":
                    missing_rep_val = value
                else:
                    na_rep_val = value

        self.missing_rep_val = missing_rep_val
        self.na_rep_val = na_rep_val

    def __len__(self):
        length = 0

        for i, variant_payload in enumerate(self.variants_payloads):
            length += consts.VARIANTS_PAYLOAD_SIZES_LEN + len(variant_payload)

            row_payload = self.variants_row_ids_payloads[i]
            if row_payload is not None:
                length += consts.ROW_IDS_SIZE_LEN + len(row_payload)

            col_payload = self.variants_col_ids_payloads[i]
            if col_payload is not None:
                length += consts.COL_IDS_SIZE_LEN + len(col_payload)

        if self.variants_amax_payload is not None:
            length += (
                consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN
                + len(self.variants_amax_payload)
            )

        if self.phase_payload is not None:
            length += consts.PHASE_PAYLOAD_SIZE_LEN + len(self.phase_payload)

            if self.phase_row_ids_payload is not None:
                length += consts.ROW_IDS_SIZE_LEN + len(self.phase_row_ids_payload)

            if self.phase_col_ids_payload is not None:
                length += consts.COL_IDS_SIZE_LEN + len(self.phase_col_ids_payload)

        if self.missing_rep_val is not None:
            length += consts.MISSING_REP_VAL_LEN
        if self.na_rep_val is not None:
            length += consts.NA_REP_VAL_LEN

        return length

    def stat(self):
        stats = {
            "NumAlleleBinMat": len(self.variants_payloads),
            "AlleleBinMat": 0,
            "AlleleRowIds": 0,
            "AlleleColIds": 0,
            "AMax": 0,
            "PhaseBinMat": 0,
            "PhaseRowIds": 0,
            "PhaseColIds": 0,
            "MissingRepVal": 0,
            "NaRepVal": 0,
        }

        for i, variant_payload in enumerate(self.variants_payloads):
            stats["AlleleBinMat"] += (
                consts.VARIANTS_PAYLOAD_SIZES_LEN + len(variant_payload)
            )

            row_payload = self.variants_row_ids_payloads[i]
            if row_payload is not None:
                stats["AlleleRowIds"] += consts.ROW_IDS_SIZE_LEN + len(row_payload)

            col_payload = self.variants_col_ids_payloads[i]
            if col_payload is not None:
                stats["AlleleColIds"] += consts.COL_IDS_SIZE_LEN + len(col_payload)

        if self.variants_amax_payload is not None:
            stats["AMax"] = (
                consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN
                + len(self.variants_amax_payload)
            )

        if self.phase_payload is not None:
            stats["PhaseBinMat"] = consts.PHASE_PAYLOAD_SIZE_LEN + len(
                self.phase_payload
            )
            if self.phase_row_ids_payload is not None:
                stats["PhaseRowIds"] = consts.ROW_IDS_SIZE_LEN + len(
                    self.phase_row_ids_payload
                )
            if self.phase_col_ids_payload is not None:
                stats["PhaseColIds"] = consts.COL_IDS_SIZE_LEN + len(
                    self.phase_col_ids_payload
                )

        if self.missing_rep_val is not None:
            stats["MissingRepVal"] = consts.MISSING_REP_VAL_LEN
        if self.na_rep_val is not None:
            stats["NaRepVal"] = consts.NA_REP_VAL_LEN

        return stats

    def to_barray(self):
        payload = bytearray()

        for i, variant_payload in enumerate(self.variants_payloads):
            variant_bytes = _payload_bytes(variant_payload)
            payload += utils.int2bstr(
                len(variant_bytes), consts.VARIANTS_PAYLOAD_SIZES_LEN
            )
            payload += variant_bytes

            row_payload = self.variants_row_ids_payloads[i]
            if row_payload is not None:
                row_bytes = _payload_bytes(row_payload)
                payload += utils.int2bstr(len(row_bytes), consts.ROW_IDS_SIZE_LEN)
                payload += row_bytes

            col_payload = self.variants_col_ids_payloads[i]
            if col_payload is not None:
                col_bytes = _payload_bytes(col_payload)
                payload += utils.int2bstr(len(col_bytes), consts.COL_IDS_SIZE_LEN)
                payload += col_bytes

        if self.variants_amax_payload is not None:
            amax_bytes = _payload_bytes(self.variants_amax_payload)
            payload += utils.int2bstr(
                len(amax_bytes), consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN
            )
            payload += amax_bytes

        if self.phase_payload is not None:
            phase_bytes = _payload_bytes(self.phase_payload)
            payload += utils.int2bstr(len(phase_bytes), consts.PHASE_PAYLOAD_SIZE_LEN)
            payload += phase_bytes

            if self.phase_row_ids_payload is not None:
                row_bytes = _payload_bytes(self.phase_row_ids_payload)
                payload += utils.int2bstr(len(row_bytes), consts.ROW_IDS_SIZE_LEN)
                payload += row_bytes

            if self.phase_col_ids_payload is not None:
                col_bytes = _payload_bytes(self.phase_col_ids_payload)
                payload += utils.int2bstr(len(col_bytes), consts.COL_IDS_SIZE_LEN)
                payload += col_bytes

        if self.missing_rep_val is not None:
            payload += utils.int2bstr(
                self.missing_rep_val, consts.MISSING_REP_VAL_LEN
            )
        if self.na_rep_val is not None:
            payload += utils.int2bstr(self.na_rep_val, consts.NA_REP_VAL_LEN)

        return payload

    def to_bytes(self):
        return bytes(self.to_barray())

    @classmethod
    def from_bitstream(cls, reader, param_set, block_payload_size):
        if block_payload_size < 0:
            raise ValueError("block payload size must be non-negative")

        start_pos = reader.tell()
        block_end = start_pos + block_payload_size

        def require_bytes(count, label):
            if count < 0 or reader.tell() + count > block_end:
                raise ValueError("truncated {} in genotype payload".format(label))

        def read_size(size_len, label):
            require_bytes(size_len, label + " size")
            size = reader.read_bytes(size_len, ret_int=True)
            require_bytes(size, label)
            return size

        def payload_region(size_len, label):
            size = read_size(size_len, label)
            # RandomAccessHandler keeps the payload lazy, but seek() itself can
            # legally move beyond physical EOF. Validate the region before
            # creating the lazy view so truncated files cannot masquerade as
            # complete blocks.
            reader.require_available(size)
            region = RandomAccessHandler(reader, reader.tell(), size)
            reader.seek(size, 1)
            return region

        variants_payloads = []
        variants_row_ids_payloads = []
        variants_col_ids_payloads = []

        for i in range(param_set.num_variants_flags):
            variants_payloads.append(
                payload_region(
                    consts.VARIANTS_PAYLOAD_SIZES_LEN,
                    "variant payload {}".format(i),
                )
            )

            if param_set.sort_variants_row_flags[i]:
                row_payload = payload_region(
                    consts.ROW_IDS_SIZE_LEN,
                    "variant row permutation {}".format(i),
                )
            else:
                row_payload = None
            variants_row_ids_payloads.append(row_payload)

            if param_set.sort_variants_col_flags[i]:
                col_payload = payload_region(
                    consts.COL_IDS_SIZE_LEN,
                    "variant column permutation {}".format(i),
                )
            else:
                col_payload = None
            variants_col_ids_payloads.append(col_payload)

        variants_amax_payload = None
        if param_set.binarization_id == consts.BinarizationID.ROW_BIN_SPLIT:
            variants_amax_payload = payload_region(
                consts.VARIANTS_AMAX_PAYLOAD_SIZE_LEN,
                "AMax payload",
            )

        phase_payload = None
        phase_row_ids_payload = None
        phase_col_ids_payload = None

        if param_set.encode_phase_data:
            phase_payload = payload_region(
                consts.PHASE_PAYLOAD_SIZE_LEN,
                "phase payload",
            )

            if param_set.sort_phases_row_flag:
                phase_row_ids_payload = payload_region(
                    consts.ROW_IDS_SIZE_LEN,
                    "phase row permutation",
                )

            if param_set.sort_phases_col_flag:
                phase_col_ids_payload = payload_region(
                    consts.COL_IDS_SIZE_LEN,
                    "phase column permutation",
                )

        if param_set.any_missing_flag:
            require_bytes(consts.MISSING_REP_VAL_LEN, "missing representation")
            missing_rep_val = reader.read_bytes(
                consts.MISSING_REP_VAL_LEN, ret_int=True
            )
        else:
            missing_rep_val = None

        if param_set.not_available_flag:
            require_bytes(consts.NA_REP_VAL_LEN, "not-available representation")
            na_rep_val = reader.read_bytes(consts.NA_REP_VAL_LEN, ret_int=True)
        else:
            na_rep_val = None

        consumed = reader.tell() - start_pos
        if consumed != block_payload_size:
            raise ValueError(
                "genotype payload length mismatch: expected {}, consumed {}".format(
                    block_payload_size, consumed
                )
            )

        payload = cls(
            param_set,
            variants_payloads,
            variants_row_ids_payloads,
            variants_col_ids_payloads,
            variants_amax_payload=variants_amax_payload,
            phase_payload=phase_payload,
            phase_row_ids_payload=phase_row_ids_payload,
            phase_col_ids_payload=phase_col_ids_payload,
            missing_rep_val=missing_rep_val,
            na_rep_val=na_rep_val,
        )

        if len(payload) != block_payload_size:
            raise ValueError("genotype payload metadata does not match serialized size")

        return payload
