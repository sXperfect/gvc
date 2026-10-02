import numpy as np

from gvc import cdebinarize, cquery, debinarize
from gvc.binarization import bin_row_bin_split
from gvc.data_structures import RowColIds, crc_id


def _query_reference(query, ploidy):
    if len(query) == 0:
        return np.empty(0, dtype=np.uint32)
    return np.concatenate(
        [
            np.arange(sample * ploidy, sample * ploidy + ploidy, dtype=np.uint32)
            for sample in query
        ]
    )


def test_cquery_expands_sample_ids_for_ploidy():
    query = np.array([1, 3], dtype=np.uint32)
    expanded = cquery.cget_col_ids(query, 2)
    np.testing.assert_array_equal(
        expanded,
        np.array([2, 3, 6, 7], dtype=np.uint32),
    )


def test_cquery_matches_reference_across_ploidies_and_query_sizes():
    rng = np.random.default_rng(3817)

    for ploidy in (1, 2, 3, 4, 8):
        for size in (0, 1, 2, 7, 31):
            query = rng.integers(0, 1000, size=size, dtype=np.uint32)
            native = cquery.cget_col_ids(query, ploidy)
            reference = _query_reference(query, ploidy)
            np.testing.assert_array_equal(native, reference)


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


def test_cython_row_split_decoder_matches_reference_for_random_matrices():
    rng = np.random.default_rng(3818)

    for nrows in (1, 2, 7, 17):
        for ncols in (1, 3, 8):
            for max_value in (0, 1, 3, 7, 15, 255):
                if max_value == 0:
                    matrix = np.zeros((nrows, ncols), dtype=np.uint8)
                else:
                    matrix = rng.integers(
                        0,
                        max_value + 1,
                        size=(nrows, ncols),
                        dtype=np.uint8,
                    )

                encoded, bit_lengths = bin_row_bin_split(matrix)
                bit_lengths_u8 = bit_lengths.astype(np.uint8)

                native = cdebinarize.debin_rc_bin_split(encoded, bit_lengths_u8)
                reference = debinarize.debin_rc_bin_split(encoded, bit_lengths_u8)

                np.testing.assert_array_equal(native, reference)
                np.testing.assert_array_equal(native, matrix)


def test_cython_permutation_decoder_matches_python_encoder():
    permutation = np.array([2, 0, 3, 1], dtype=np.uint16)
    payload = RowColIds(permutation).to_bytes()
    encoded = np.frombuffer(payload, dtype=np.uint8)

    restored = crc_id.decode_rowcolids(encoded, len(permutation))

    np.testing.assert_array_equal(restored, permutation)


def test_cython_permutation_decoder_matches_python_reference_randomized():
    rng = np.random.default_rng(3819)

    for size in (2, 3, 4, 7, 8, 15, 16, 31, 64, 255):
        for _ in range(8):
            permutation = rng.permutation(size).astype(np.uint16)
            payload = RowColIds(permutation).to_bytes()
            encoded = np.frombuffer(payload, dtype=np.uint8)

            native = crc_id.decode_rowcolids(encoded, size)
            reference = RowColIds.from_bytes(payload, size).ids

            np.testing.assert_array_equal(native, reference)
            np.testing.assert_array_equal(native, permutation)



def _expected_uniform_phase_text(allele_mat, ploidy, phase_value):
    separator = "|" if phase_value == 0 else "/"
    rows = []
    n_samples = allele_mat.shape[1] // ploidy
    for row in allele_mat:
        genotypes = []
        for sample in range(n_samples):
            values = row[sample * ploidy:(sample + 1) * ploidy]
            genotypes.append(
                separator.join("." if int(value) == -1 else str(int(value)) for value in values)
            )
        rows.append("\t".join(genotypes))
    return "\n".join(rows) + "\n"


def test_native_uniform_phase_text_matches_gvc_phase_convention():
    rng = np.random.default_rng(3820)

    for ploidy in (1, 2, 3, 4):
        for phase_value in (0, 1):
            matrix = rng.integers(
                0,
                5,
                size=(5, 4 * ploidy),
                dtype=np.int8,
            )
            matrix[0, 0] = -1

            expected = _expected_uniform_phase_text(matrix, ploidy, phase_value)
            python_text = debinarize.recon_gt_mat_with_phase_val(
                matrix.copy(), phase_value, ploidy
            )
            native_text = cdebinarize.recon_gt_mat_with_phase_val(
                matrix.copy(), phase_value, ploidy
            ).decode("ascii")

            assert python_text == expected
            assert native_text == expected


def test_python_mixed_phase_text_uses_zero_for_pipe_one_for_slash():
    allele_mat = np.array([[0, 1, 1, 0]], dtype=np.int8)
    phase_mat = np.array([[0, 1]], dtype=np.uint8)
    assert debinarize.recon_gt_mat_with_phase_mat(
        allele_mat.copy(), phase_mat, 2
    ) == "0|1\t1/0\n"


def test_native_row_split_rejects_inconsistent_bit_lengths():
    encoded = np.zeros((2, 3), dtype=bool)
    with np.testing.assert_raises_regex(
        ValueError, "does not match encoded row count"
    ):
        cdebinarize.debin_rc_bin_split(
            encoded, np.array([1], dtype=np.uint8)
        )


def test_cquery_rejects_nonpositive_ploidy():
    query = np.array([1, 3], dtype=np.uint32)
    with np.testing.assert_raises_regex(ValueError, "greater than zero"):
        cquery.cget_col_ids(query, 0)



def test_ctypes_libgvc_permutation_decoder_matches_python_reference(monkeypatch):
    from pathlib import Path

    from gvc import libds

    library = Path(__file__).resolve().parents[1] / "tmp" / "libgvc-build" / "libds.so"
    if not library.is_file():
        import pytest
        pytest.skip("standalone libgvc has not been built")

    monkeypatch.setenv("GVC_LIBDS_PATH", str(library))
    monkeypatch.setattr(libds, "_LIBDS", None)

    rng = np.random.default_rng(3821)
    for size in (1, 2, 3, 4, 7, 8, 15, 16, 31, 255):
        permutation = rng.permutation(size).astype(np.uint16)
        payload = RowColIds(permutation).to_bytes()
        restored = libds.decode_rowcolids(payload, size)
        np.testing.assert_array_equal(restored, permutation)


def test_ctypes_libgvc_rejects_truncated_and_invalid_payload(monkeypatch):
    from pathlib import Path

    import pytest

    from gvc import libds

    library = Path(__file__).resolve().parents[1] / "tmp" / "libgvc-build" / "libds.so"
    if not library.is_file():
        pytest.skip("standalone libgvc has not been built")

    monkeypatch.setenv("GVC_LIBDS_PATH", str(library))
    monkeypatch.setattr(libds, "_LIBDS", None)

    permutation = np.array([2, 0, 3, 1], dtype=np.uint16)
    payload = RowColIds(permutation).to_bytes()
    with pytest.raises(ValueError, match="truncated"):
        libds.decode_rowcolids(payload[:-1], len(permutation))

    # For three entries, two bits are used per id; value 3 is outside 0..2.
    with pytest.raises(ValueError, match="valid range"):
        libds.decode_rowcolids(bytes([0b11000000]), 3)



def test_ctypes_libgvc_rejects_duplicate_and_trailing_payload(monkeypatch):
    from pathlib import Path

    import pytest

    from gvc import libds

    library = (
        Path(__file__).resolve().parents[1]
        / "tmp"
        / "libgvc-build"
        / "libds.so"
    )
    if not library.is_file():
        pytest.skip("standalone libgvc has not been built")

    monkeypatch.setenv("GVC_LIBDS_PATH", str(library))
    monkeypatch.setattr(libds, "_LIBDS", None)

    with pytest.raises(ValueError, match="duplicate"):
        libds.decode_rowcolids(bytes([0b00000100]), 3)

    permutation = np.array([2, 0, 3, 1], dtype=np.uint16)
    payload = RowColIds(permutation).to_bytes()
    with pytest.raises(ValueError, match="trailing"):
        libds.decode_rowcolids(payload + b"\x00", len(permutation))



def test_python_genotype_text_supports_multi_digit_alleles():
    from gvc import debinarize

    allele_mat = np.array([[10, 11, 2, 12]], dtype=np.int8)
    assert debinarize.recon_gt_mat_with_phase_val(
        allele_mat.copy(), 0, 2
    ) == "10|11\t2|12\n"


def test_python_genotype_text_supports_multi_digit_mixed_phase():
    from gvc import debinarize

    allele_mat = np.array([[10, 11, 12, 13]], dtype=np.int8)
    phase_mat = np.array([[0, 1]], dtype=np.uint8)
    assert debinarize.recon_gt_mat_with_phase_mat(
        allele_mat.copy(), phase_mat, 2
    ) == "10|11\t12/13\n"


def test_python_genotype_text_rejects_inconsistent_phase_shape():
    from gvc import debinarize

    allele_mat = np.array([[0, 1, 1, 0]], dtype=np.int8)
    with pytest.raises(ValueError, match="phasing tensor shape"):
        debinarize.simd_tensor_to_txt(
            allele_mat.reshape(1, 2, 2),
            np.zeros((1, 1, 1), dtype=np.uint8),
        )



def test_native_row_split_rejects_bit_length_above_uint8_width():
    encoded = np.zeros((9, 1), dtype=bool)
    with np.testing.assert_raises_regex(
        ValueError, "must not exceed 8"
    ):
        cdebinarize.debin_rc_bin_split(
            encoded, np.array([9], dtype=np.uint8)
        )



def test_cquery_rejects_uint32_expansion_overflow():
    query = np.array([np.iinfo(np.uint32).max], dtype=np.uint32)
    with np.testing.assert_raises_regex(
        ValueError, "exceeds uint32 range"
    ):
        cquery.cget_col_ids(query, 2)



def test_native_helpers_reject_ploidy_above_format_limit():
    query = np.array([0], dtype=np.uint32)
    with np.testing.assert_raises_regex(
        ValueError, "serialized GVC limit"
    ):
        cquery.cget_col_ids(query, 257)

    matrix = np.zeros((1, 257), dtype=np.int8)
    with np.testing.assert_raises_regex(
        ValueError, "serialized GVC limit"
    ):
        cdebinarize.recon_gt_mat_with_phase_val(matrix, 0, 257)
