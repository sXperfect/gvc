import numpy as np

from gvc import cdebinarize, cquery
from gvc import debinarize
from gvc.binarization import bin_row_bin_split
from gvc.data_structures import RowColIds
from gvc.data_structures import crc_id


def test_cquery_expands_sample_ids_for_ploidy():
    query = np.array([1, 3], dtype=np.uint32)
    expanded = cquery.cget_col_ids(query, 2)
    np.testing.assert_array_equal(
        expanded,
        np.array([2, 3, 6, 7], dtype=np.uint32),
    )


def test_cython_row_split_decoder_matches_python_reference():
    matrix = np.array(
        [
            [0, 1, 2, 3],
            [7, 0, 4, 1],
            [0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    encoded, bit_lengths = bin_row_bin_split(matrix)
    bit_lengths_u8 = bit_lengths.astype(np.uint8)

    native = cdebinarize.debin_rc_bin_split(encoded, bit_lengths_u8)
    reference = debinarize.debin_rc_bin_split(encoded, bit_lengths_u8)

    np.testing.assert_array_equal(native, matrix)
    np.testing.assert_array_equal(native, reference)


def test_cython_permutation_decoder_matches_python_encoder():
    permutation = np.array([2, 0, 3, 1], dtype=np.uint16)
    payload = RowColIds(permutation).to_bytes()
    encoded = np.frombuffer(payload, dtype=np.uint8)

    restored = crc_id.decode_rowcolids(encoded, len(permutation))

    np.testing.assert_array_equal(restored, permutation)
