from pathlib import Path

import numpy as np
import pytest

from gvc import binarization, reader


pytest.importorskip("cyvcf2")

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_diploid.vcf"
HAPLOID_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_haploid.vcf"
MIXED_PLOIDY_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_mixed_ploidy.vcf"


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



def test_haploid_vcf_reader_has_zero_width_phase_matrix():
    blocks = list(
        reader.vcf_genotypes_reader(str(HAPLOID_FIXTURE), None, block_size=2)
    )
    assert len(blocks) == 2

    alleles, phases, ploidy, missing_rep, na_rep = blocks[0]
    assert ploidy == 1
    assert phases.shape == (2, 0)
    assert missing_rep is None
    assert na_rep is None
    np.testing.assert_array_equal(
        alleles,
        np.array([[0, 1], [2, 0]], dtype=np.uint8),
    )

    tail_alleles, tail_phases, ploidy, missing_rep, na_rep = blocks[1]
    assert ploidy == 1
    assert tail_phases.shape == (1, 0)
    assert missing_rep == 2
    assert na_rep is None
    np.testing.assert_array_equal(
        tail_alleles,
        np.array([[2, 1]], dtype=np.uint8),
    )



def test_exact_block_boundary_finalizes_metadata(tmp_path):
    output = tmp_path / "exact.gvc"
    blocks = list(
        reader.vcf_genotypes_reader(
            str(FIXTURE),
            str(output),
            block_size=3,
        )
    )
    assert len(blocks) == 1

    metadata = Path(str(output) + ".metadata")
    np.testing.assert_array_equal(
        np.load(metadata / "main.npy"),
        np.array([[100, 300]], dtype=np.uint64),
    )
    np.testing.assert_array_equal(
        np.load(metadata / "0.npy"),
        np.array([100, 200, 300], dtype=np.uint64),
    )


def test_reader_splits_blocks_when_ploidy_changes():
    blocks = list(
        reader.vcf_genotypes_reader(
            str(MIXED_PLOIDY_FIXTURE),
            None,
            block_size=4,
        )
    )
    assert [block[2] for block in blocks] == [1, 2]

    haploid_alleles, haploid_phases, _, _, _ = blocks[0]
    np.testing.assert_array_equal(
        haploid_alleles,
        np.array([[0, 1], [2, 0]], dtype=np.uint8),
    )
    assert haploid_phases.shape == (2, 0)

    diploid_alleles, diploid_phases, _, _, _ = blocks[1]
    np.testing.assert_array_equal(
        diploid_alleles,
        np.array([[0, 1, 1, 1], [2, 1, 0, 2]], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        diploid_phases,
        np.array([[False, True], [True, False]], dtype=bool),
    )



def test_metadata_end_writes_empty_root_index(tmp_path):
    class FakeVCF:
        raw_header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A"
        )
        samples = ["SAMPLE_A"]

    metadata = tmp_path / "empty.gvc.metadata"
    handler = reader.MetaHandler(FakeVCF(), str(metadata), block_size=4)
    handler.init()
    handler.end()

    root = np.load(metadata / "main.npy", allow_pickle=False)
    assert root.shape == (0, 2)
    assert root.dtype == np.uint64



def test_vcf_reader_rejects_allele_index_that_overflows_working_dtype(monkeypatch):
    class FakeGenotype:
        def array(self):
            return np.array([[128, 0, 1]], dtype=np.int16)

    class FakeVariant:
        ploidy = 2
        POS = 100
        genotype = FakeGenotype()

    class FakeVCF:
        raw_header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A"
        )
        samples = ["SAMPLE_A"]

        def __init__(self, *args, **kwargs):
            self._variants = [FakeVariant()]

        def __iter__(self):
            return iter(self._variants)

        def close(self):
            pass

    monkeypatch.setattr(reader, "_vcf_class", lambda: FakeVCF)

    with pytest.raises(ValueError, match="outside the supported range"):
        list(reader.vcf_genotypes_reader("synthetic.vcf", None, block_size=1))



def test_metadata_samples_use_authoritative_vcf_sample_list(tmp_path):
    class FakeVCF:
        raw_header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWRONG"
        )
        samples = ["SAMPLE_A", "SAMPLE_B"]

    metadata = tmp_path / "samples.gvc.metadata"
    handler = reader.MetaHandler(FakeVCF(), str(metadata), block_size=4)
    handler.init()

    np.testing.assert_array_equal(
        np.load(metadata / "samples.npy", allow_pickle=False),
        np.array(["SAMPLE_A", "SAMPLE_B"]),
    )
