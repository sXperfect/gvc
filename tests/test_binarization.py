import numpy as np
import pytest

from gvc import binarization


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_bit_plane_roundtrip(axis):
    matrix = np.array(
        [
            [0, 1, 2, 3],
            [3, 2, 1, 0],
        ],
        dtype=np.uint8,
    )
    encoded, bit_depth = binarization.bin_bit_plane(matrix, axis=axis)
    decoded = binarization.debin_bit_plane(encoded, bit_depth, axis)
    np.testing.assert_array_equal(decoded, matrix)


def test_bit_plane_all_zero_matrix_has_one_plane():
    matrix = np.zeros((3, 4), dtype=np.uint8)
    encoded, bit_depth = binarization.bin_bit_plane(matrix, axis=2)
    assert bit_depth == 1
    assert len(encoded) == 1
    np.testing.assert_array_equal(
        binarization.debin_bit_plane(encoded, bit_depth, 2),
        matrix,
    )


def test_row_bin_split_roundtrip():
    matrix = np.array(
        [
            [0, 1, 2, 3],
            [7, 0, 4, 1],
            [0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    encoded, bit_lengths = binarization.bin_row_bin_split(matrix)
    decoded = binarization.debin_row_bin_split([encoded], bit_lengths)
    np.testing.assert_array_equal(decoded, matrix)


def test_adaptive_missing_values_roundtrip():
    source = np.array([[0, -1, -2, 3], [2, 1, -1, 0]], dtype=np.int8)
    encoded, missing_value, na_value = binarization.adaptive_max_value(source)
    assert missing_value is not None
    assert na_value is not None
    decoded = binarization.undo_adaptive_max_value(
        encoded.copy(), missing_value, na_value
    )
    np.testing.assert_array_equal(decoded, source)


def test_split_genotype_matrix_handles_phase_and_missing_values():
    lines = [
        "0|1\t1/1\n",
        "./.\t2|0\n",
    ]
    allele_matrix, phase_matrix, ploidy = binarization.split_genotype_matrix(lines)

    assert ploidy == 2
    np.testing.assert_array_equal(
        allele_matrix,
        np.array(
            [
                [0, 1, 1, 1],
                [-1, -1, 2, 0],
            ],
            dtype=np.int8,
        ),
    )
    np.testing.assert_array_equal(
        phase_matrix,
        np.array(
            [
                [0, 1],
                [1, 0],
            ],
            dtype=bool,
        ),
    )
