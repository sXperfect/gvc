import numpy as np

from gvc.data_structures import RowColIds, VectorAMax


def test_row_col_ids_roundtrip():
    permutation = np.array([2, 0, 3, 1], dtype=np.uint16)
    payload = RowColIds(permutation).to_bitio().to_bytes(align=True)
    restored = RowColIds.from_bytes(payload, len(permutation)).ids
    np.testing.assert_array_equal(restored, permutation)


def test_singleton_row_col_ids_roundtrip():
    permutation = np.array([0], dtype=np.uint16)
    payload = RowColIds(permutation).to_bitio().to_bytes(align=True)
    restored = RowColIds.from_bytes(payload, 1).ids
    np.testing.assert_array_equal(restored, permutation)


def test_amax_roundtrip():
    values = np.array([1, 2, 3, 7, 1], dtype=np.uint16)
    payload = VectorAMax(values).to_bitio().to_bytes(align=True)
    restored = VectorAMax.from_bytes(payload).vector
    np.testing.assert_array_equal(restored, values)


def test_empty_row_col_ids_roundtrip():
    permutation = np.array([], dtype=np.uint16)
    payload = RowColIds(permutation).to_bitio().to_bytes(align=True)
    restored = RowColIds.from_bytes(payload, 0).ids
    np.testing.assert_array_equal(restored, permutation)



def test_amax_rejects_impossible_entry_count_before_allocation():
    # Header claims 2^32-1 entries but provides no entry flags.
    payload = (2**32 - 1).to_bytes(4, "big") + b"\x01"

    with np.testing.assert_raises_regex(
        ValueError, "entry count exceeds available data"
    ):
        VectorAMax.from_bytes(payload)



def test_amax_rejects_value_above_uint8_bit_length_limit():
    with np.testing.assert_raises_regex(
        ValueError, "must not exceed 8"
    ):
        VectorAMax(np.array([9], dtype=np.uint16))


def test_amax_rejects_serialized_entry_width_above_format_limit():
    payload = (1).to_bytes(4, "big") + b"\x04" + b"\x00"

    with np.testing.assert_raises_regex(
        ValueError, "entry width exceeds"
    ):
        VectorAMax.from_bytes(payload)



def test_row_col_decoder_rejects_entry_count_above_uint16_range():
    with np.testing.assert_raises_regex(
        ValueError, "outside the uint16 permutation range"
    ):
        RowColIds.from_bytes(b"", np.iinfo(np.uint16).max + 1)
