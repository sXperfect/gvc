import json
from pathlib import Path

import pytest

from scripts import run_release_validation


def test_release_validation_rejects_threshold_without_baseline(tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"x")

    with pytest.raises(SystemExit) as exc_info:
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--max-regression-percent",
                "5",
            ]
        )
    assert exc_info.value.code == 2


def test_release_validation_rejects_duplicate_workers(tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"x")

    with pytest.raises(SystemExit) as exc_info:
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--workers",
                "1",
                "1",
            ]
        )
    assert exc_info.value.code == 2


def test_release_validation_writes_machine_readable_evidence(monkeypatch, tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    benchmark = tmp_path / "benchmark.json"
    evidence = tmp_path / "evidence.json"

    commands = []

    def fake_run(command, env=None):
        commands.append(list(command))
        if "benchmarks/run_release.py" in command:
            benchmark.write_text(
                json.dumps({"schema_version": 1, "configurations": []}),
                encoding="utf-8",
            )

    monkeypatch.setattr(run_release_validation, "_run", fake_run)

    status = run_release_validation.main(
        [
            str(fixture),
            "--benchmark-output",
            str(benchmark),
            "--evidence-output",
            str(evidence),
            "--workers",
            "0",
            "2",
            "--repetitions",
            "1",
            "--historical-max-blocks",
            "1",
        ]
    )

    assert status == 0
    payload = json.loads(evidence.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 2
    assert payload["status"] == "pass"
    assert payload["commit"] == run_release_validation._git_commit()
    assert isinstance(payload["git_dirty"], bool)
    assert payload["package_version"]
    assert payload["host"]["python"]
    assert payload["host"]["executable"]
    assert payload["fixture"]["bytes"] == len(b"fixture")
    assert len(payload["fixture"]["sha256"]) == 64
    assert payload["fixture"]["git_blob"]
    assert payload["configuration"] == {
        "historical_max_blocks": 1,
        "include_sorting": False,
        "workers": [0, 2],
        "repetitions": 1,
        "block_size": 2048,
        "start_method": None,
        "max_regression_percent": None,
    }
    assert [step["name"] for step in payload["steps"]] == [
        "static_release_readiness",
        "historical_fixture",
        "release_benchmark",
    ]
    assert any("scripts/check_release_readiness.py" in command for command in commands)
    assert any("scripts/verify_historical.py" in command for command in commands)
    assert any("benchmarks/run_release.py" in command for command in commands)


def test_rc_validation_requires_full_historical_fixture(tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"x")

    with pytest.raises(SystemExit) as exc_info:
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--rc",
                "--historical-max-blocks",
                "1",
                "--include-sorting",
            ]
        )
    assert exc_info.value.code == 2


def test_rc_validation_requires_sorting_coverage(tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"x")

    with pytest.raises(SystemExit) as exc_info:
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--rc",
                "--historical-max-blocks",
                "0",
            ]
        )
    assert exc_info.value.code == 2


def test_release_validation_persists_failure_evidence(monkeypatch, tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    benchmark = tmp_path / "benchmark.json"
    evidence = tmp_path / "evidence.json"

    def fake_run(command, env=None):
        if "scripts/verify_historical.py" in command:
            raise RuntimeError("synthetic historical failure")

    monkeypatch.setattr(run_release_validation, "_run", fake_run)

    with pytest.raises(RuntimeError, match="synthetic historical failure"):
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(benchmark),
                "--evidence-output",
                str(evidence),
                "--workers",
                "0",
                "--repetitions",
                "1",
                "--historical-max-blocks",
                "1",
            ]
        )

    payload = json.loads(evidence.read_text(encoding="utf-8"))
    assert payload["status"] == "fail"
    assert payload["steps"][-1]["name"] == "historical_fixture"
    assert payload["steps"][-1]["status"] == "fail"
    assert "synthetic historical failure" in payload["steps"][-1]["error"]
    assert payload["completed_unix"] >= payload["created_unix"]


def test_rc_artifact_failure_persists_evidence(monkeypatch, tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    evidence = tmp_path / "evidence.json"

    monkeypatch.setattr(run_release_validation, "_git_dirty", lambda: False)
    monkeypatch.setattr(run_release_validation, "_package_version", lambda: "1.0.1rc1")

    def fake_run(command, env=None):
        if "-m" in command and "build" in command:
            raise RuntimeError("synthetic artifact build failure")

    monkeypatch.setattr(run_release_validation, "_run", fake_run)

    with pytest.raises(RuntimeError, match="synthetic artifact build failure"):
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--evidence-output",
                str(evidence),
                "--rc",
                "--historical-max-blocks",
                "0",
                "--include-sorting",
            ]
        )

    payload = json.loads(evidence.read_text(encoding="utf-8"))
    assert payload["status"] == "fail"
    assert payload["steps"][-1]["name"] == "release_artifacts"
    assert payload["steps"][-1]["status"] == "fail"
    assert payload["steps"][-1]["version"] == "1.0.1rc1"
    assert "synthetic artifact build failure" in payload["steps"][-1]["error"]


def test_rc_baseline_requires_explicit_regression_budget(tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    baseline = tmp_path / "baseline.json"
    baseline.write_text('{"schema_version": 1, "configurations": []}', encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--baseline",
                str(baseline),
                "--rc",
                "--historical-max-blocks",
                "0",
                "--include-sorting",
            ]
        )
    assert exc_info.value.code == 2


def test_rc_dirty_worktree_persists_failure_evidence(monkeypatch, tmp_path):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    evidence = tmp_path / "evidence.json"

    monkeypatch.setattr(run_release_validation, "_git_dirty", lambda: True)
    monkeypatch.setattr(run_release_validation, "_package_version", lambda: "1.0.1rc1")

    with pytest.raises(RuntimeError, match="clean git worktree"):
        run_release_validation.main(
            [
                str(fixture),
                "--benchmark-output",
                str(tmp_path / "benchmark.json"),
                "--evidence-output",
                str(evidence),
                "--rc",
                "--historical-max-blocks",
                "0",
                "--include-sorting",
            ]
        )

    payload = json.loads(evidence.read_text(encoding="utf-8"))
    assert payload["status"] == "fail"
    assert payload["steps"][-1]["name"] == "clean_worktree"
    assert payload["steps"][-1]["status"] == "fail"
    assert "clean git worktree" in payload["steps"][-1]["error"]



def test_release_validation_requires_complete_metric_comparison(
    monkeypatch,
    tmp_path,
):
    fixture = tmp_path / "fixture.vcf.gz"
    fixture.write_bytes(b"fixture")
    benchmark = tmp_path / "benchmark.json"
    baseline = tmp_path / "baseline.json"
    evidence = tmp_path / "evidence.json"
    baseline.write_text(
        json.dumps({"schema_version": 1, "configurations": []}),
        encoding="utf-8",
    )

    commands = []

    def fake_run(command, env=None):
        commands.append(list(command))
        if "benchmarks/run_release.py" in command:
            benchmark.write_text(
                json.dumps({"schema_version": 1, "configurations": []}),
                encoding="utf-8",
            )
        if "benchmarks/compare_release.py" in command:
            output_index = command.index("--output") + 1
            Path(command[output_index]).write_text(
                json.dumps({}),
                encoding="utf-8",
            )

    monkeypatch.setattr(run_release_validation, "_run", fake_run)

    assert run_release_validation.main(
        [
            str(fixture),
            "--benchmark-output",
            str(benchmark),
            "--evidence-output",
            str(evidence),
            "--baseline",
            str(baseline),
            "--max-regression-percent",
            "5",
            "--workers",
            "0",
            "--repetitions",
            "1",
            "--historical-max-blocks",
            "1",
        ]
    ) == 0

    compare_commands = [
        command
        for command in commands
        if "benchmarks/compare_release.py" in command
    ]
    assert len(compare_commands) == 1
    assert "--require-complete-metrics" in compare_commands[0]
