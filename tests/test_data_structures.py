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
