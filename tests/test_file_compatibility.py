from io import BytesIO
from pathlib import Path

import pytest

from gvc.bitstream import BitstreamReader
from gvc.data_structures import AccessUnit, ParameterSet
from gvc.data_structures.consts import DataUnitType
from gvc.data_structures.data_unit import DataUnitHeader


FIXTURE = Path(__file__).parent / "fixtures" / "v1_structural.gvc"


def _parse_structural_v1(data):
    reader = BitstreamReader(BytesIO(data))

    first_type = reader.read_bytes(1, ret_int=True)
    if first_type != DataUnitType.PARAMETER_SET:
        raise ValueError("expected parameter set as first data unit")
    parameter_header = DataUnitHeader.from_bitstream(first_type, reader)
    parameter_set = ParameterSet.from_bitstream(reader, parameter_header)

    second_type = reader.read_bytes(1, ret_int=True)
    if second_type != DataUnitType.ACCESS_UNIT:
        raise ValueError("expected access unit after parameter set")
    access_unit = AccessUnit.from_bitstream(
        reader,
        {parameter_set.parameter_set_id: parameter_set},
    )

    if reader.tell() != len(data):
        raise ValueError("trailing bytes after access unit")

    return parameter_set, access_unit


def test_complete_v1_structural_fixture_is_stable():
    data = FIXTURE.read_bytes()
    assert data.hex() == (
        "000000000a0000401800"
        "0100000016000000070001000000000600000002aa55"
    )

    parameter_set, access_unit = _parse_structural_v1(data)

    assert parameter_set.parameter_set_id == 0
    assert parameter_set.p == 2
    assert parameter_set.num_bin_mat == 1
    assert parameter_set.concat_axis == 2

    assert access_unit.header.access_unit_id == 7
    assert access_unit.header.parameter_set_id == 0
    assert access_unit.num_blocks == 1
    assert (
        access_unit.blocks[0].block_payload.variants_payloads[0].read()
        == b"\xaa\x55"
    )


@pytest.mark.parametrize("cut", range(32))
def test_complete_v1_fixture_rejects_every_truncation(cut):
    data = FIXTURE.read_bytes()[:cut]
    with pytest.raises((EOFError, ValueError)):
        _parse_structural_v1(data)


def test_complete_v1_fixture_rejects_unknown_access_unit_type():
    data = bytearray(FIXTURE.read_bytes())
    data[10] = 0x7F
    with pytest.raises(ValueError, match="expected access unit"):
        _parse_structural_v1(data)


def test_complete_v1_fixture_rejects_unknown_block_content_id():
    data = bytearray(FIXTURE.read_bytes())
    data[21] = 0x01
    with pytest.raises(ValueError, match="unsupported block content id"):
        _parse_structural_v1(data)


def test_complete_v1_fixture_rejects_declared_access_unit_length_mismatch():
    data = bytearray(FIXTURE.read_bytes())
    data[14] += 1
    with pytest.raises((EOFError, ValueError)):
        _parse_structural_v1(data)


def test_complete_v1_fixture_rejects_trailing_bytes():
    with pytest.raises(ValueError, match="trailing bytes"):
        _parse_structural_v1(FIXTURE.read_bytes() + b"\x00")
