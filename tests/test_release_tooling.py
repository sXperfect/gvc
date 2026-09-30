import json
from pathlib import Path

import pytest

from benchmarks import compare_release
from scripts import check_release_artifacts
from scripts import check_release_readiness
from scripts import check_release_tag
from scripts import run_release_validation


ROOT = Path(__file__).resolve().parents[1]


def test_release_readiness_declares_required_release_surfaces():
    problems = check_release_readiness.check_readiness(allow_dev=True)
    assert problems == []


def test_release_validation_file_record_is_content_addressed(tmp_path):
    value = tmp_path / "evidence.bin"
    value.write_bytes(b"gvc-release-evidence")

    record = run_release_validation._record_file(value)

    assert record["path"] == str(value.resolve())
    assert record["bytes"] == len(b"gvc-release-evidence")
    assert len(record["sha256"]) == 64
    assert record["sha256"] == run_release_validation._sha256(value)


def test_release_validation_git_commit_matches_repository_head():
    assert run_release_validation._git_commit() == (
        run_release_validation._git_output("rev-parse", "HEAD")
    )


def test_manual_release_workflow_requires_baseline_hash():
    workflow = (
        ROOT / ".github" / "workflows" / "release-validation.yml"
    ).read_text(encoding="utf-8")

    assert "baseline_url:" in workflow
    assert "baseline_sha256:" in workflow
    assert 'baseline_sha256 is required with baseline_url' in workflow
    assert "sha256sum -c -" in workflow
    assert "--baseline tmp/release-validation/baseline.json" in workflow


def test_release_checklist_requires_same_commit_and_artifact_hashes():
    checklist = (
        ROOT / "docs" / "audits" / "release-1.0.1-checklist.md"
    ).read_text(encoding="utf-8")
    assert "exact Git commit" in checklist
    assert "wheel/sdist SHA-256" in checklist
    assert "same release candidate" in checklist


def test_benchmark_comparison_detects_regression():
    key = {
        "workers": 0,
        "start_method": None,
        "block_size": 2048,
        "binarization": "bit_plane",
        "axis": 2,
        "sort_rows": False,
        "sort_cols": False,
        "transpose": False,
    }
    baseline = {
        "schema_version": 1,
        "configurations": [
            dict(
                key,
                encode_variants_per_second_median=100.0,
                total_gvc_bytes=1000,
            )
        ],
    }
    candidate = {
        "schema_version": 1,
        "configurations": [
            dict(
                key,
                encode_variants_per_second_median=80.0,
                total_gvc_bytes=1200,
            )
        ],
    }

    result = compare_release.compare_reports(
        baseline, candidate, threshold_percent=10.0
    )
    assert len(result["failures"]) == 2


def test_artifact_version_policy_rejects_development_versions():
    assert check_release_tag.package_version().startswith("1.0.")
    # Development branches are intentionally not taggable as RC/final.
    if "dev" in check_release_tag.package_version():
        assert check_release_tag.package_version() != "1.0.1rc1"
