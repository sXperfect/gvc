#!/usr/bin/env python3
"""Create reproducible GVC release benchmark reports."""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
WORKER = ROOT / "benchmarks" / "_release_worker.py"


def _run_configuration(args, workers):
    command = [
        sys.executable,
        str(WORKER),
        args.fixture,
        "--workers",
        str(workers),
        "--block-size",
        str(args.block_size),
        "--binarization",
        args.binarization,
        "--axis",
        str(args.axis),
    ]
    if args.sort_rows:
        command.append("--sort-rows")
    if args.sort_cols:
        command.append("--sort-cols")
    if args.transpose:
        command.append("--transpose")
    if args.start_method:
        command += ["--start-method", args.start_method]
    if args.stall_timeout is not None:
        command += ["--stall-timeout", str(args.stall_timeout)]

    samples = []
    for _ in range(args.repetitions):
        completed = subprocess.run(
            command,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "benchmark worker failed with status {}:\n{}".format(
                    completed.returncode,
                    completed.stderr,
                )
            )
        samples.append(json.loads(completed.stdout.strip()))

    metric_names = (
        "encode_seconds",
        "decode_seconds",
        "random_access_seconds",
        "encode_variants_per_second",
        "decode_variants_per_second",
        "peak_rss_kib_self",
        "peak_rss_kib_children",
    )
    summary = dict(samples[-1])
    summary["repetitions"] = len(samples)
    summary["samples"] = samples
    for metric in metric_names:
        values = [sample[metric] for sample in samples if sample[metric] is not None]
        if values:
            summary[metric + "_median"] = statistics.median(values)
            summary[metric + "_min"] = min(values)
            summary[metric + "_max"] = max(values)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run isolated GVC release benchmark configurations."
    )
    parser.add_argument("fixture")
    parser.add_argument("--workers", type=int, nargs="+", default=[0, 1, 2, 4])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument("--binarization", choices=["bit_plane", "row_bin_split"], default="bit_plane")
    parser.add_argument("--axis", type=int, choices=[0, 1, 2], default=2)
    parser.add_argument("--sort-rows", action="store_true")
    parser.add_argument("--sort-cols", action="store_true")
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--start-method", choices=["fork", "spawn", "forkserver"])
    parser.add_argument("--stall-timeout", type=float)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    if args.repetitions <= 0:
        parser.error("--repetitions must be positive")
    if any(value < 0 for value in args.workers):
        parser.error("--workers values must be non-negative")
    if len(set(args.workers)) != len(args.workers):
        parser.error("--workers values must be unique")
    if not os.path.isfile(args.fixture):
        parser.error("fixture does not exist: {}".format(args.fixture))

    report = {
        "schema_version": 1,
        "created_unix": time.time(),
        "fixture": os.path.abspath(args.fixture),
        "fixture_bytes": os.path.getsize(args.fixture),
        "host": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cpu_count": os.cpu_count(),
        },
        "configurations": [
            _run_configuration(args, workers)
            for workers in args.workers
        ],
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
