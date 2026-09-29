from io import BytesIO
from pathlib import Path

import pytest

from gvc.bitstream import BitstreamReader
from gvc.data_structures import AccessUnit, ParameterSet
from gvc.data_structures.consts import BinarizationID, CodecID, DataUnitType
from gvc.data_structures.data_unit import DataUnitHeader


FIXTURES = Path(__file__).parent / "fixtures"

CASES = {
    "v1_row_bin_missing.gvc": {
        "hex": (
            "000000000901c06620"
            "010000002d000000090101000000001d"
            "0000000311223300000001900000000140"
            "000000060000000302d00405"
        ),
        "parameter_set_id": 1,
        "access_unit_id": 9,
        "binarization": BinarizationID.ROW_BIN_SPLIT,
    },
    "v1_mixed_phase.gvc": {
        "hex": (
            "000000000b0200401a0c00"
            "01000000240000000a02010000000014"
            "000000015500000001aa00000001800000000100"
        ),
        "parameter_set_id": 2,
        "access_unit_id": 10,
        "binarization": BinarizationID.BIT_PLANE,
    },
    "v1_multi_plane_sorted.gvc": {
        "hex": (
            "000000000b030040290680"
            "01000000250000000b03010000000015"
            "000000016100000001800000000262630000000140"
        ),
        "parameter_set_id": 3,
        "access_unit_id": 11,
        "binarization": BinarizationID.BIT_PLANE,
    },
}


def _parse(data):
    reader = BitstreamReader(BytesIO(data))

    data_unit_type = reader.read_bytes(1, ret_int=True)
    if data_unit_type != DataUnitType.PARAMETER_SET:
        raise ValueError("expected parameter set")
    header = DataUnitHeader.from_bitstream(data_unit_type, reader)
    parameter_set = ParameterSet.from_bitstream(reader, header)

    data_unit_type = reader.read_bytes(1, ret_int=True)
    if data_unit_type != DataUnitType.ACCESS_UNIT:
        raise ValueError("expected access unit")
    access_unit = AccessUnit.from_bitstream(
        reader,
        {parameter_set.parameter_set_id: parameter_set},
    )

    if reader.tell() != len(data):
        raise ValueError("trailing bytes after access unit")
    return parameter_set, access_unit


@pytest.mark.parametrize("name", sorted(CASES))
def test_v1_structural_corpus_is_byte_exact_and_decodable(name):
    case = CASES[name]
    data = (FIXTURES / name).read_bytes()

    assert data.hex() == case["hex"]
    parameter_set, access_unit = _parse(data)

    assert parameter_set.parameter_set_id == case["parameter_set_id"]
    assert parameter_set.binarization_id == case["binarization"]
    assert access_unit.header.access_unit_id == case["access_unit_id"]
    assert access_unit.header.parameter_set_id == case["parameter_set_id"]
    assert access_unit.num_blocks == 1


def test_row_bin_missing_fixture_preserves_optional_payload_order():
    _, access_unit = _parse((FIXTURES / "v1_row_bin_missing.gvc").read_bytes())
    payload = access_unit.blocks[0].block_payload

    assert payload.variants_payloads[0].read() == b"\x11\x22\x33"
    assert payload.variants_row_ids_payloads[0].read() == b"\x90"
    assert payload.variants_col_ids_payloads[0].read() == b"\x40"
    assert payload.variants_amax_payload.read() == bytes.fromhex(
        "0000000302d0"
    )
    assert payload.missing_rep_val == 4
    assert payload.na_rep_val == 5


def test_mixed_phase_fixture_preserves_phase_payload_order():
    parameter_set, access_unit = _parse(
        (FIXTURES / "v1_mixed_phase.gvc").read_bytes()
    )
    payload = access_unit.blocks[0].block_payload

    assert parameter_set.encode_phase_data is True
    assert parameter_set.sort_phases_row_flag is True
    assert parameter_set.sort_phases_col_flag is True
    assert payload.variants_payloads[0].read() == b"\x55"
    assert payload.phase_payload.read() == b"\xaa"
    assert payload.phase_row_ids_payload.read() == b"\x80"
    assert payload.phase_col_ids_payload.read() == b"\x00"


def test_multi_plane_fixture_preserves_per_plane_flags_and_payloads():
    parameter_set, access_unit = _parse(
        (FIXTURES / "v1_multi_plane_sorted.gvc").read_bytes()
    )
    payload = access_unit.blocks[0].block_payload

    assert parameter_set.num_bin_mat == 2
    assert parameter_set.concat_axis == 2
    assert parameter_set.sort_variants_row_flags == [True, False]
    assert parameter_set.sort_variants_col_flags == [False, True]
    assert parameter_set.transpose_variants_mat_flags == [False, True]
    assert parameter_set.variants_coder_ids == [CodecID.JBIG1, CodecID.GABAC]
    assert payload.variants_payloads[0].read() == b"a"
    assert payload.variants_row_ids_payloads[0].read() == b"\x80"
    assert payload.variants_payloads[1].read() == b"bc"
    assert payload.variants_col_ids_payloads[1].read() == b"\x40"


@pytest.mark.parametrize("name", sorted(CASES))
def test_v1_structural_corpus_rejects_every_truncation(name):
    data = (FIXTURES / name).read_bytes()
    for cut in range(len(data)):
        with pytest.raises((EOFError, ValueError, TypeError)):
            _parse(data[:cut])
