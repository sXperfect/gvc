#!/usr/bin/env python3
"""Orchestrate controlled GVC 1.0 release-candidate validation."""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def _run(command, env=None):
    print("$ " + " ".join(str(value) for value in command), flush=True)
    completed = subprocess.run(
        [str(value) for value in command],
        cwd=str(ROOT),
        env=env,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            "command failed with status {}: {}".format(
                completed.returncode,
                " ".join(str(value) for value in command),
            )
        )


def _package_version():
    namespace = {}
    version_file = ROOT / "gvc" / "_version.py"
    exec(version_file.read_text(encoding="utf-8"), namespace)
    return namespace["__version__"]


def _artifact(path):
    value = Path(path)
    if not value.is_file():
        raise FileNotFoundError("required file does not exist: {}".format(value))
    return value.resolve()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run controlled GVC 1.0 release-candidate validation."
    )
    parser.add_argument("fixture", help="historical LUH VCF fixture")
    parser.add_argument(
        "--benchmark-output",
        required=True,
        help="path for the candidate benchmark JSON report",
    )
    parser.add_argument(
        "--baseline",
        help="optional retained benchmark JSON to compare against",
    )
    parser.add_argument(
        "--max-regression-percent",
        type=float,
        default=None,
        help="optional benchmark regression budget",
    )
    parser.add_argument(
        "--workers",
        type=int,
        nargs="+",
        default=[0, 1, 2, 4],
    )
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument(
        "--historical-max-blocks",
        type=int,
        default=0,
        help="0 processes the full historical fixture",
    )
    parser.add_argument(
        "--include-sorting",
        action="store_true",
        help="include all row/column sorting combinations in historical checks",
    )
    parser.add_argument(
        "--start-method",
        choices=["fork", "spawn", "forkserver"],
    )
    parser.add_argument(
        "--evidence-output",
        help="optional JSON summary of completed validation steps",
    )
    parser.add_argument(
        "--rc",
        action="store_true",
        help="require an RC/final package version rather than .dev",
    )
    args = parser.parse_args(argv)

    fixture = _artifact(args.fixture)
    if args.baseline:
        baseline = _artifact(args.baseline)
    else:
        baseline = None

    if args.repetitions <= 0:
        parser.error("--repetitions must be positive")
    if args.block_size <= 0:
        parser.error("--block-size must be positive")
    if args.historical_max_blocks < 0:
        parser.error("--historical-max-blocks must be non-negative")
    if any(worker < 0 for worker in args.workers):
        parser.error("--workers values must be non-negative")
    if len(set(args.workers)) != len(args.workers):
        parser.error("--workers values must be unique")
    if (
        args.max_regression_percent is not None
        and args.max_regression_percent < 0
    ):
        parser.error("--max-regression-percent must be non-negative")
    if args.max_regression_percent is not None and baseline is None:
        parser.error("--max-regression-percent requires --baseline")

    evidence = {
        "schema_version": 1,
        "created_unix": time.time(),
        "fixture": str(fixture),
        "fixture_bytes": fixture.stat().st_size,
        "rc_mode": bool(args.rc),
        "steps": [],
    }

    readiness = [sys.executable, "scripts/check_release_readiness.py"]
    if args.rc:
        readiness.append("--rc")
    _run(readiness)
    evidence["steps"].append({"name": "static_release_readiness", "status": "pass"})

    if args.rc:
        dist_dir = ROOT / "tmp" / "release-validation" / "dist"
        if dist_dir.exists():
            import shutil
            shutil.rmtree(str(dist_dir))
        dist_dir.mkdir(parents=True, exist_ok=True)
        _run([
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--sdist",
            "--outdir",
            str(dist_dir),
        ])
        version = _package_version()
        wheels = sorted(dist_dir.glob("gvc-*.whl"))
        sdists = sorted(dist_dir.glob("gvc-*.tar.gz"))
        if len(wheels) != 1 or len(sdists) != 1:
            raise RuntimeError("expected exactly one wheel and one sdist")
        artifact_evidence = (
            ROOT / "tmp" / "release-validation" / "artifacts.json"
        )
        _run([
            sys.executable,
            "scripts/check_release_artifacts.py",
            "--wheel",
            str(wheels[0]),
            "--sdist",
            str(sdists[0]),
            "--version",
            version,
            "--output",
            str(artifact_evidence),
        ])
        evidence["steps"].append(
            {
                "name": "release_artifacts",
                "status": "pass",
                "version": version,
                "evidence": str(artifact_evidence),
            }
        )

    historical = [
        sys.executable,
        "scripts/verify_historical.py",
        str(fixture),
        "--block-size",
        str(args.block_size),
        "--max-blocks",
        str(args.historical_max_blocks),
    ]
    if args.include_sorting:
        historical.append("--include-sorting")
    _run(historical)
    evidence["steps"].append(
        {
            "name": "historical_fixture",
            "status": "pass",
            "max_blocks": args.historical_max_blocks,
            "include_sorting": bool(args.include_sorting),
        }
    )

    benchmark_output = Path(args.benchmark_output).resolve()
    benchmark = [
        sys.executable,
        "benchmarks/run_release.py",
        str(fixture),
        "--workers",
    ]
    benchmark.extend(str(worker) for worker in args.workers)
    benchmark += [
        "--repetitions",
        str(args.repetitions),
        "--block-size",
        str(args.block_size),
        "--output",
        str(benchmark_output),
    ]
    if args.start_method:
        benchmark += ["--start-method", args.start_method]
    _run(benchmark)
    evidence["steps"].append(
        {
            "name": "release_benchmark",
            "status": "pass",
            "report": str(benchmark_output),
            "workers": args.workers,
            "repetitions": args.repetitions,
            "start_method": args.start_method,
        }
    )

    if baseline is not None:
        compare = [
            sys.executable,
            "benchmarks/compare_release.py",
            str(baseline),
            str(benchmark_output),
            "--require-same-configurations",
        ]
        if args.max_regression_percent is not None:
            compare += [
                "--max-regression-percent",
                str(args.max_regression_percent),
            ]
        _run(compare)
        evidence["steps"].append(
            {
                "name": "benchmark_comparison",
                "status": "pass",
                "baseline": str(baseline),
                "max_regression_percent": args.max_regression_percent,
            }
        )

    evidence["status"] = "pass"
    if args.evidence_output:
        output = Path(args.evidence_output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(evidence, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print("release evidence: {}".format(output))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
