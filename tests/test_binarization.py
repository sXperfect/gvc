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



def test_split_genotype_matrix_handles_haploid_input():
    allele_matrix, phase_matrix, ploidy = binarization.split_genotype_matrix(
        ["0\t1\n", ".\t2\n"]
    )
    assert ploidy == 1
    assert phase_matrix is None
    np.testing.assert_array_equal(
        allele_matrix,
        np.array([[0, 1], [-1, 2]], dtype=np.int8),
    )


def test_adaptive_max_value_rejects_unsigned_input():
    with pytest.raises(TypeError, match="signed integer"):
        binarization.adaptive_max_value(
            np.array([[0, 1]], dtype=np.uint8)
        )


def test_binarizers_reject_invalid_shapes_and_values():
    with pytest.raises(ValueError, match="two-dimensional"):
        binarization.bin_bit_plane(np.zeros((2, 2, 1), dtype=np.uint8), axis=2)
    with pytest.raises(ValueError, match="non-negative"):
        binarization.bin_bit_plane(np.array([[0, -1]], dtype=np.int8), axis=2)
    with pytest.raises(ValueError, match="axis"):
        binarization.bin_bit_plane(np.array([[0, 1]], dtype=np.uint8), axis=7)

    with pytest.raises(ValueError, match="two-dimensional"):
        binarization.bin_row_bin_split(np.zeros((2, 2, 1), dtype=np.uint8))


def test_debin_bit_plane_rejects_plane_count_mismatch():
    plane = np.zeros((2, 2), dtype=bool)
    with pytest.raises(ValueError, match="count"):
        binarization.debin_bit_plane([plane], bit_depth=2, axis=2)
    with pytest.raises(ValueError, match="one matrix"):
        binarization.debin_bit_plane([plane, plane], bit_depth=2, axis=0)



def test_adaptive_max_value_accepts_reader_tensor():
    source = np.array(
        [
            [[0, -1], [1, 1]],
            [[-2, 2], [0, -1]],
        ],
        dtype=np.int8,
    )
    encoded, missing_value, na_value = binarization.adaptive_max_value(source)
    assert encoded.shape == source.shape
    restored = binarization.undo_adaptive_max_value(
        encoded.copy(), missing_value, na_value
    )
    np.testing.assert_array_equal(restored, source)



def test_adaptive_missing_roundtrip_at_int8_upper_boundary():
    source = np.array([[127, -1, -2, 0]], dtype=np.int8)

    encoded, missing_value, na_value = binarization.adaptive_max_value(source)

    assert int(missing_value) == 128
    assert int(na_value) == 129
    restored = binarization.undo_adaptive_max_value(
        encoded, missing_value, na_value
    )
    np.testing.assert_array_equal(restored, source)



def test_matrix_to_tensor_rejects_invalid_grouping():
    matrix = np.zeros((2, 5), dtype=np.uint8)

    with pytest.raises(ValueError, match="positive"):
        binarization.matrix_to_tensor(matrix, 0)

    with pytest.raises(ValueError, match="divisible"):
        binarization.matrix_to_tensor(matrix, 2)



def test_split_genotype_matrix_rejects_allele_above_int8_range():
    with pytest.raises(ValueError, match="supported range"):
        binarization.split_genotype_matrix(["128|0\n"])
