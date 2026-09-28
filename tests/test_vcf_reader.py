from pathlib import Path

import numpy as np
import pytest

from gvc import binarization, reader


pytest.importorskip("cyvcf2")

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_diploid.vcf"


def test_tiny_vcf_reader_preserves_alleles_and_gvc_phase_convention():
    blocks = list(reader.vcf_genotypes_reader(str(FIXTURE), None, block_size=2))
    assert len(blocks) == 2

    alleles, phases, ploidy, missing_rep, na_rep = blocks[0]
    assert ploidy == 2
    assert missing_rep is None
    assert na_rep is None
    np.testing.assert_array_equal(
        alleles,
        np.array(
            [
                [0, 1, 1, 1],
                [2, 1, 0, 2],
            ],
            dtype=np.uint8,
        ),
    )
    # GVC convention: 0 means "|" (phased), 1 means "/" (unphased).
    np.testing.assert_array_equal(
        phases,
        np.array(
            [
                [False, True],
                [True, False],
            ],
            dtype=bool,
        ),
    )

    tail_alleles, tail_phases, ploidy, missing_rep, na_rep = blocks[1]
    assert ploidy == 2
    assert missing_rep == 2
    assert na_rep is None
    np.testing.assert_array_equal(
        tail_alleles,
        np.array([[2, 2, 1, 0]], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        tail_phases,
        np.array([[True, False]], dtype=bool),
    )


def test_tiny_vcf_reader_matches_text_parser_for_complete_records():
    lines = ["0|1\t1/1\n", "2/1\t0|2\n"]
    expected_alleles, expected_phases, expected_ploidy = (
        binarization.split_genotype_matrix(lines)
    )

    alleles, phases, ploidy, _, _ = next(
        reader.vcf_genotypes_reader(str(FIXTURE), None, block_size=2)
    )

    assert ploidy == expected_ploidy
    np.testing.assert_array_equal(alleles, expected_alleles)
    np.testing.assert_array_equal(phases, expected_phases)
