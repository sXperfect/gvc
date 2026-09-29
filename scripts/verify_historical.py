#!/usr/bin/env python3
"""Run historical LUH fixture compatibility checks with the production codec."""

import argparse
import itertools
import os
import sys

import numpy as np

from gvc.codec import jbigkit
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.decoder import decode_encoded_variants
from gvc.encoder import run_core
from gvc.reader import vcf_genotypes_reader


def _cases(include_sorting):
    if include_sorting:
        sort_modes = list(itertools.product((False, True), repeat=2))
    else:
        sort_modes = [(False, False)]

    for binarization_id, axis in (
        (BinarizationID.BIT_PLANE, 2),
        (BinarizationID.ROW_BIN_SPLIT, 0),
    ):
        for sort_rows, sort_cols in sort_modes:
            yield binarization_id, axis, sort_rows, sort_cols


def verify_case(
    fixture,
    block_size,
    max_blocks,
    binarization_id,
    axis,
    sort_rows,
    sort_cols,
):
    ps_params = [
        binarization_id,
        CodecID.JBIG1,
        axis,
        sort_rows,
        sort_cols,
        False,
    ]
    tsp_params = ["ham", "nn", 0]

    seen = 0
    for raw_block in vcf_genotypes_reader(fixture, None, block_size):
        original_alleles = raw_block[0].copy()
        original_phase = raw_block[1].copy()

        block, parameter_set = run_core(raw_block, ps_params, tsp_params)
        restored_alleles, restored_phase = decode_encoded_variants(
            parameter_set,
            block.block_payload,
            ret_gt=False,
        )

        np.testing.assert_array_equal(restored_alleles, original_alleles)
        if parameter_set.p == 1:
            if restored_phase not in (0, False):
                raise AssertionError("unexpected haploid phase value")
        elif parameter_set.encode_phase_data:
            np.testing.assert_array_equal(restored_phase, original_phase)
        elif original_phase.size:
            if bool(restored_phase) != bool(original_phase.flat[0]):
                raise AssertionError("constant phase value changed")

        seen += 1
        if max_blocks is not None and seen >= max_blocks:
            break

    if seen == 0:
        raise RuntimeError("fixture did not yield any genotype blocks")
    if max_blocks is not None and seen < max_blocks:
        raise RuntimeError(
            "fixture yielded {} blocks, expected at least {}".format(
                seen, max_blocks
            )
        )
    return seen


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Verify GVC against the original LUH VCF fixture."
    )
    parser.add_argument("fixture", help="path to test_block01.vcf.gz or equivalent")
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument(
        "--max-blocks",
        type=int,
        default=1,
        help="number of blocks per case; 0 means all blocks",
    )
    parser.add_argument(
        "--include-sorting",
        action="store_true",
        help="also exercise all row/column sorting combinations",
    )
    parser.add_argument(
        "--encoder",
        help="override pbmtojbg85 path (or use GVC_JBIG_ENCODER)",
    )
    parser.add_argument(
        "--decoder",
        help="override jbgtopbm85 path (or use GVC_JBIG_DECODER)",
    )
    args = parser.parse_args(argv)

    if args.block_size <= 0:
        parser.error("--block-size must be positive")
    if args.max_blocks < 0:
        parser.error("--max-blocks must be non-negative")
    if not os.path.isfile(args.fixture):
        parser.error("fixture does not exist: {}".format(args.fixture))

    if args.encoder:
        os.environ[jbigkit.ENCODER_ENV] = args.encoder
    if args.decoder:
        os.environ[jbigkit.DECODER_ENV] = args.decoder

    if not jbigkit.available():
        parser.error(
            "JBIG-KIT executables are unavailable; install jbigkit-bin or "
            "set GVC_JBIG_ENCODER/GVC_JBIG_DECODER"
        )

    max_blocks = None if args.max_blocks == 0 else args.max_blocks
    total = 0
    for binarization_id, axis, sort_rows, sort_cols in _cases(
        args.include_sorting
    ):
        count = verify_case(
            args.fixture,
            args.block_size,
            max_blocks,
            binarization_id,
            axis,
            sort_rows,
            sort_cols,
        )
        total += count
        print(
            "PASS binarization={} axis={} sort_rows={} sort_cols={} blocks={}".format(
                binarization_id.name,
                axis,
                sort_rows,
                sort_cols,
                count,
            )
        )

    print("historical compatibility checks passed: {} block-cases".format(total))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
