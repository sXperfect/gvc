#!/usr/bin/env python3
"""Run one isolated GVC release-benchmark configuration."""

import argparse
import json
import os
import platform
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from gvc.codec import jbigkit
from gvc.decoder import Decoder
from gvc.encoder import Encoder


def _metadata_bytes(metadata_dir):
    total = 0
    for path in Path(metadata_dir).rglob("*"):
        if path.is_file():
            total += path.stat().st_size
    return total


def _variant_count(metadata_dir):
    metadata = Path(metadata_dir)
    main = np.load(metadata / "main.npy", allow_pickle=False)
    total = 0
    for block_id in range(main.shape[0]):
        positions = np.load(metadata / ("{}.npy".format(block_id)), allow_pickle=False)
        total += int(positions.shape[0])
    return total


def _first_sample(metadata_dir):
    samples = np.load(Path(metadata_dir) / "samples.npy", allow_pickle=False)
    if samples.size == 0:
        return None
    return str(samples[0])


def _rss_kib(who):
    value = resource.getrusage(who).ru_maxrss
    if sys.platform == "darwin":
        return int(value // 1024)
    return int(value)


def run_once(args):
    jbigkit.install()

    with tempfile.TemporaryDirectory(prefix="gvc-benchmark-") as directory:
        root = Path(directory)
        encoded = root / "benchmark.gvc"
        decoded = root / "decoded.txt"
        selected = root / "selected.txt"

        encoder = Encoder(
            args.fixture,
            str(encoded),
            binarization_name=args.binarization,
            axis=args.axis,
            sort_rows=args.sort_rows,
            sort_cols=args.sort_cols,
            transpose=args.transpose,
            block_size=args.block_size,
            codec_name="jbig",
            preset_mode=0,
            num_threads=args.workers,
            multiprocessing_start_method=args.start_method,
            multiprocessing_stall_timeout=args.stall_timeout,
            multiprocessing_initializer=(
                jbigkit.install if args.workers > 0 else None
            ),
        )

        start = time.perf_counter()
        encoder.run()
        encode_seconds = time.perf_counter() - start

        decoder = Decoder(str(encoded), str(decoded))
        try:
            start = time.perf_counter()
            decoder.decode()
            decode_seconds = time.perf_counter() - start
        finally:
            if decoder._out_f is not None:
                decoder._out_f.close()
            decoder._f.close()

        metadata_dir = Path(str(encoded) + ".metadata")
        variants = _variant_count(metadata_dir)
        sample = _first_sample(metadata_dir)

        random_access_seconds = None
        if sample is not None:
            ra_decoder = Decoder(str(encoded), str(selected))
            try:
                start = time.perf_counter()
                ra_decoder.random_access(None, sample)
                random_access_seconds = time.perf_counter() - start
            finally:
                if ra_decoder._out_f is not None:
                    ra_decoder._out_f.close()
                ra_decoder._f.close()

        gvc_stream_bytes = encoded.stat().st_size
        metadata_bytes = _metadata_bytes(metadata_dir)
        total_gvc_bytes = gvc_stream_bytes + metadata_bytes
        decoded_gt_bytes = decoded.stat().st_size

        return {
            "workers": args.workers,
            "start_method": args.start_method,
            "block_size": args.block_size,
            "binarization": args.binarization,
            "axis": args.axis,
            "sort_rows": args.sort_rows,
            "sort_cols": args.sort_cols,
            "transpose": args.transpose,
            "variants": variants,
            "encode_seconds": encode_seconds,
            "decode_seconds": decode_seconds,
            "random_access_seconds": random_access_seconds,
            "encode_variants_per_second": (
                variants / encode_seconds if encode_seconds > 0 else None
            ),
            "decode_variants_per_second": (
                variants / decode_seconds if decode_seconds > 0 else None
            ),
            "gvc_stream_bytes": gvc_stream_bytes,
            "metadata_bytes": metadata_bytes,
            "total_gvc_bytes": total_gvc_bytes,
            "decoded_gt_bytes": decoded_gt_bytes,
            "decoded_gt_to_total_gvc_ratio": (
                decoded_gt_bytes / total_gvc_bytes
                if total_gvc_bytes
                else None
            ),
            "peak_rss_kib_self": _rss_kib(resource.RUSAGE_SELF),
            "peak_rss_kib_children": _rss_kib(resource.RUSAGE_CHILDREN),
            "python": platform.python_version(),
            "platform": platform.platform(),
        }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture")
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument("--binarization", choices=["bit_plane", "row_bin_split"], default="bit_plane")
    parser.add_argument("--axis", type=int, choices=[0, 1, 2], default=2)
    parser.add_argument("--sort-rows", action="store_true")
    parser.add_argument("--sort-cols", action="store_true")
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--start-method", choices=["fork", "spawn", "forkserver"])
    parser.add_argument("--stall-timeout", type=float)
    args = parser.parse_args(argv)

    if args.workers < 0:
        parser.error("--workers must be non-negative")
    if args.block_size <= 0:
        parser.error("--block-size must be positive")
    if not os.path.isfile(args.fixture):
        parser.error("fixture does not exist: {}".format(args.fixture))
    if args.binarization == "row_bin_split" and args.axis != 0:
        parser.error("row_bin_split benchmark requires --axis 0")

    print(json.dumps(run_once(args), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
