import numpy as np
import pytest

from gvc.bitstream import BitIO, BitstreamWriter
from gvc.data_structures import AccessUnit, Block
from gvc.data_structures.amax import VectorAMax
from gvc.data_structures.access_unit import AccessUnitHeader
from gvc.data_structures.block import BlockHeader
from gvc.data_structures.data_unit import DataUnitHeader
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.data_structures.param_set import ParameterSet
from gvc.data_structures.rc_id import RowColIds


def _parameter_set(**overrides):
    values = dict(
        parameter_set_id=0,
        any_missing_flag=False,
        not_available_flag=False,
        p=2,
        binarization_id=BinarizationID.BIT_PLANE,
        num_bin_mat=1,
        concat_axis=2,
        sort_variants_row_flags=[False],
        sort_variants_col_flags=[False],
        transpose_variants_mat_flags=[False],
        variants_coder_ids=[CodecID.JBIG1],
        encode_phase_data=False,
        phase_value=False,
    )
    values.update(overrides)
    return ParameterSet(**values)


@pytest.mark.parametrize("p", [1, 2, 255, 256])
def test_ploidy_serialized_boundary_accepts_valid_values(p):
    parameter_set = _parameter_set(p=p)
    assert parameter_set.p == p


@pytest.mark.parametrize("p", [0, 257])
def test_ploidy_serialized_boundary_rejects_invalid_values(p):
    with pytest.raises(ValueError, match="serialized range"):
        _parameter_set(p=p)


@pytest.mark.parametrize("num_bin_mat", [1, 2, 254, 255])
def test_bitplane_count_serialized_boundary_accepts_valid_values(num_bin_mat):
    flags = [False] * num_bin_mat
    parameter_set = _parameter_set(
        num_bin_mat=num_bin_mat,
        concat_axis=2,
        sort_variants_row_flags=flags,
        sort_variants_col_flags=flags,
        transpose_variants_mat_flags=flags,
        variants_coder_ids=[CodecID.JBIG1] * num_bin_mat,
    )
    assert parameter_set.num_bin_mat == num_bin_mat


def test_bitplane_count_rejects_unserializable_value():
    with pytest.raises(ValueError, match="num_bin_mat"):
        _parameter_set(num_bin_mat=256)


@pytest.mark.parametrize("axis", [-1, 3, 4])
def test_concat_axis_rejects_out_of_range_values(axis):
    with pytest.raises(ValueError, match="concat_axis"):
        _parameter_set(concat_axis=axis)


def test_permutation_rejects_id_outside_domain():
    with pytest.raises(ValueError, match="valid range"):
        RowColIds(np.array([0, 4, 1, 2], dtype=np.uint16)).to_bytes()


def test_permutation_decoder_rejects_truncated_payload():
    with pytest.raises(ValueError, match="truncated"):
        RowColIds.from_bytes(b"", 4)


def test_amax_rejects_zero_and_non_vector_inputs():
    with pytest.raises(ValueError, match="non-zero"):
        VectorAMax(np.array([1, 0, 2], dtype=np.uint16))
    with pytest.raises(ValueError, match="one-dimensional"):
        VectorAMax(np.ones((2, 2), dtype=np.uint16))


def test_bitio_rejects_values_that_do_not_fit():
    bits = BitIO()
    with pytest.raises(RuntimeError, match="exceeds"):
        bits.write(4, 2)


def test_block_header_rejects_negative_payload_size():
    with pytest.raises(ValueError, match="non-negative"):
        BlockHeader(0, -1)


def test_block_requires_genotype_payload_type():
    with pytest.raises(TypeError):
        Block.from_encoded_variant(object())



def test_serialized_header_ranges_are_explicit():
    with pytest.raises(ValueError, match="data-unit type"):
        DataUnitHeader(256, 0)
    with pytest.raises(ValueError, match="data-unit length"):
        DataUnitHeader(0, 1 << 32)

    with pytest.raises(ValueError, match="access_unit_id"):
        AccessUnitHeader(6, 1 << 32, 0, 0)
    with pytest.raises(ValueError, match="parameter_set_id"):
        AccessUnitHeader(6, 0, 256, 0)
    with pytest.raises(ValueError, match="num_blocks"):
        AccessUnitHeader(6, 0, 0, 256)

    with pytest.raises(ValueError, match="payload size"):
        BlockHeader(0, 1 << 32)


def test_access_unit_rejects_non_block_payload():
    header = AccessUnitHeader(6, 0, 0, 1)
    with pytest.raises(TypeError, match="Block"):
        AccessUnit(header, [object()])


def test_permutation_requires_true_permutation():
    with pytest.raises(ValueError, match="unique"):
        RowColIds(np.array([0, 0, 2], dtype=np.uint16))
    with pytest.raises(ValueError, match="valid range"):
        RowColIds(np.array([1], dtype=np.uint16))


def test_permutation_rejects_trailing_bytes():
    permutation = RowColIds(np.array([2, 0, 3, 1], dtype=np.uint16))
    with pytest.raises(ValueError, match="trailing"):
        RowColIds.from_bytes(permutation.to_bytes() + b"\x00", 4)


def test_amax_accepts_empty_but_rejects_negative_values_and_trailing_bytes():
    empty = VectorAMax(np.array([], dtype=np.uint16))
    encoded = empty.to_bitio().to_bytes(align=True)
    assert VectorAMax.from_bytes(encoded).vector.size == 0

    with pytest.raises(ValueError, match="positive"):
        VectorAMax(np.array([1, -1], dtype=np.int16))

    payload = VectorAMax(np.array([1, 2, 3], dtype=np.uint16)).to_bitio().to_bytes(
        align=True
    )
    with pytest.raises(ValueError, match="trailing"):
        VectorAMax.from_bytes(payload + b"\x00")


def test_bit_writers_reject_negative_values_and_widths():
    import io

    with pytest.raises(ValueError, match="non-negative"):
        BitIO().write(-1, 1)
    with pytest.raises(ValueError, match="nbits"):
        BitIO().write(0, -1)

    writer = BitstreamWriter(io.BytesIO())
    with pytest.raises(ValueError, match="non-negative"):
        writer.write_bits(-1, 1)
    with pytest.raises(ValueError, match="nbits"):
        writer.write_bits(0, -1)
