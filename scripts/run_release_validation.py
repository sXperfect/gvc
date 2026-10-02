#!/usr/bin/env python3
"""Orchestrate controlled GVC 1.0 release-candidate validation."""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
HISTORICAL_FIXTURE_GIT_BLOB = "af45a419e46563906ac51fad0be869291316cea9"


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


def _write_evidence(path, evidence):
    if not path:
        return
    output = Path(path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("release evidence: {}".format(output))


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


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_output(*args):
    completed = subprocess.run(
        ["git"] + list(args),
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    return completed.stdout.strip()


def _git_commit():
    return _git_output("rev-parse", "HEAD")


def _git_dirty():
    return bool(_git_output("status", "--porcelain"))


def _git_blob(path):
    return _git_output("hash-object", str(path))


def _install_and_smoke_artifact(path, name):
    from scripts import ci

    status = ci.install_and_smoke_artifact(Path(path).resolve(), name)
    if status:
        raise RuntimeError(
            "{} artifact install/smoke failed with status {}".format(
                name, status
            )
        )


def _record_file(path):
    value = Path(path).resolve()
    return {
        "path": str(value),
        "bytes": value.stat().st_size,
        "sha256": _sha256(value),
    }


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
    if args.rc and baseline is not None and args.max_regression_percent is None:
        parser.error("--rc with --baseline requires --max-regression-percent")
    if args.rc and args.historical_max_blocks != 0:
        parser.error("--rc requires --historical-max-blocks 0 for the full fixture")
    if args.rc and not args.include_sorting:
        parser.error("--rc requires --include-sorting")

    evidence = {
        "schema_version": 2,
        "created_unix": time.time(),
        "commit": _git_commit(),
        "git_dirty": _git_dirty(),
        "package_version": _package_version(),
        "host": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "executable": sys.executable,
            "cpu_count": os.cpu_count(),
        },
        "fixture": dict(
            _record_file(fixture),
            git_blob=_git_blob(fixture),
        ),
        "rc_mode": bool(args.rc),
        "configuration": {
            "historical_max_blocks": args.historical_max_blocks,
            "include_sorting": bool(args.include_sorting),
            "workers": args.workers,
            "repetitions": args.repetitions,
            "block_size": args.block_size,
            "start_method": args.start_method,
            "max_regression_percent": args.max_regression_percent,
        },
        "steps": [],
    }
    if args.rc and evidence["fixture"]["git_blob"] != HISTORICAL_FIXTURE_GIT_BLOB:
        evidence["steps"].append(
            {
                "name": "historical_fixture_identity",
                "status": "fail",
                "error": (
                    "RC validation requires pinned historical fixture blob {}"
                    .format(HISTORICAL_FIXTURE_GIT_BLOB)
                ),
                "actual_git_blob": evidence["fixture"]["git_blob"],
            }
        )
        evidence["status"] = "fail"
        evidence["completed_unix"] = time.time()
        _write_evidence(args.evidence_output, evidence)
        raise RuntimeError("RC validation historical fixture identity mismatch")

    if args.rc and evidence["git_dirty"]:
        evidence["steps"].append(
            {
                "name": "clean_worktree",
                "status": "fail",
                "error": "RC validation requires a clean git worktree",
            }
        )
        evidence["status"] = "fail"
        evidence["completed_unix"] = time.time()
        _write_evidence(args.evidence_output, evidence)
        raise RuntimeError("RC validation requires a clean git worktree")

    readiness = [sys.executable, "scripts/check_release_readiness.py"]
    if args.rc:
        readiness.append("--rc")
    try:
        _run(readiness)
    except Exception as exc:
        evidence["steps"].append(
            {"name": "static_release_readiness", "status": "fail", "error": str(exc)}
        )
        evidence["status"] = "fail"
        evidence["completed_unix"] = time.time()
        _write_evidence(args.evidence_output, evidence)
        raise
    evidence["steps"].append({"name": "static_release_readiness", "status": "pass"})

    if args.rc:
        dist_dir = ROOT / "tmp" / "release-validation" / "dist"
        try:
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
            artifact_data = json.loads(
                artifact_evidence.read_text(encoding="utf-8")
            )
            _install_and_smoke_artifact(wheels[0], "release-wheel")
            _install_and_smoke_artifact(sdists[0], "release-sdist")
        except Exception as exc:
            evidence["steps"].append(
                {
                    "name": "release_artifacts",
                    "status": "fail",
                    "error": str(exc),
                    "version": _package_version(),
                }
            )
            evidence["status"] = "fail"
            evidence["completed_unix"] = time.time()
            _write_evidence(args.evidence_output, evidence)
            raise
        evidence["steps"].append(
            {
                "name": "release_artifacts",
                "status": "pass",
                "version": version,
                "evidence": _record_file(artifact_evidence),
                "artifacts": artifact_data["artifacts"],
                "install_smoke": ["wheel", "sdist"],
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
    try:
        _run(historical)
    except Exception as exc:
        evidence["steps"].append(
            {
                "name": "historical_fixture",
                "status": "fail",
                "error": str(exc),
                "max_blocks": args.historical_max_blocks,
                "include_sorting": bool(args.include_sorting),
            }
        )
        evidence["status"] = "fail"
        evidence["completed_unix"] = time.time()
        _write_evidence(args.evidence_output, evidence)
        raise
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
    try:
        _run(benchmark)
    except Exception as exc:
        evidence["steps"].append(
            {
                "name": "release_benchmark",
                "status": "fail",
                "error": str(exc),
                "workers": args.workers,
                "repetitions": args.repetitions,
                "start_method": args.start_method,
            }
        )
        evidence["status"] = "fail"
        evidence["completed_unix"] = time.time()
        _write_evidence(args.evidence_output, evidence)
        raise
    evidence["steps"].append(
        {
            "name": "release_benchmark",
            "status": "pass",
            "report": _record_file(benchmark_output),
            "workers": args.workers,
            "repetitions": args.repetitions,
            "start_method": args.start_method,
        }
    )

    if baseline is not None:
        comparison_output = (
            ROOT / "tmp" / "release-validation" / "comparison.json"
        )
        compare = [
            sys.executable,
            "benchmarks/compare_release.py",
            str(baseline),
            str(benchmark_output),
            "--require-same-configurations",
            "--require-compatible-environment",
            "--require-complete-metrics",
            "--output",
            str(comparison_output),
        ]
        if args.max_regression_percent is not None:
            compare += [
                "--max-regression-percent",
                str(args.max_regression_percent),
            ]
        try:
            _run(compare)
        except Exception as exc:
            evidence["steps"].append(
                {
                    "name": "benchmark_comparison",
                    "status": "fail",
                    "error": str(exc),
                    "baseline": _record_file(baseline),
                    "max_regression_percent": args.max_regression_percent,
                }
            )
            evidence["status"] = "fail"
            evidence["completed_unix"] = time.time()
            _write_evidence(args.evidence_output, evidence)
            raise
        evidence["steps"].append(
            {
                "name": "benchmark_comparison",
                "status": "pass",
                "baseline": _record_file(baseline),
                "comparison": _record_file(comparison_output),
                "max_regression_percent": args.max_regression_percent,
            }
        )

    evidence["status"] = "pass"
    evidence["completed_unix"] = time.time()
    _write_evidence(args.evidence_output, evidence)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
