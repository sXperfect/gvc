import json

import pytest

from benchmarks.compare_release import compare_reports


def _config(workers, encode_vps=100.0, decode_vps=200.0, size=1000):
    return {
        "workers": workers,
        "start_method": None,
        "block_size": 128,
        "binarization": "bit_plane",
        "axis": 2,
        "sort_rows": False,
        "sort_cols": False,
        "transpose": False,
        "encode_variants_per_second_median": encode_vps,
        "decode_variants_per_second_median": decode_vps,
        "encode_seconds_median": 10.0,
        "decode_seconds_median": 5.0,
        "random_access_seconds_median": 0.1,
        "peak_rss_kib_self_median": 10000,
        "peak_rss_kib_children_median": 20000,
        "total_gvc_bytes": size,
    }


def _report(configurations):
    return {"schema_version": 1, "configurations": configurations}


def test_compare_release_detects_directional_regressions():
    baseline = _report([_config(2)])
    candidate = _report([
        _config(2, encode_vps=80.0, decode_vps=220.0, size=1100)
    ])

    result = compare_reports(baseline, candidate, threshold_percent=15.0)
    failed = {entry[1]["metric"] for entry in result["failures"]}

    assert "encode_variants_per_second_median" in failed
    assert "total_gvc_bytes" not in failed
    assert "decode_variants_per_second_median" not in failed


def test_compare_release_reports_configuration_mismatch():
    result = compare_reports(
        _report([_config(0), _config(2)]),
        _report([_config(2), _config(4)]),
    )
    assert len(result["missing_configurations"]) == 1
    assert len(result["extra_configurations"]) == 1


def test_compare_release_rejects_duplicate_configuration():
    with pytest.raises(ValueError, match="duplicate benchmark configuration"):
        compare_reports(
            _report([_config(2), _config(2)]),
            _report([_config(2)]),
        )


def test_compare_release_reports_fixture_provenance_mismatch():
    baseline = _report([_config(2)])
    candidate = _report([_config(2)])
    baseline["fixture_sha256"] = "a" * 64
    candidate["fixture_sha256"] = "b" * 64
    baseline["host"] = {
        "platform": "linux",
        "machine": "x86_64",
        "processor": "cpu",
        "python": "3.8.20",
        "cpu_count": 8,
    }
    candidate["host"] = dict(baseline["host"])

    result = compare_reports(baseline, candidate)

    assert result["provenance_mismatches"] == ["fixture_sha256"]


def test_compare_release_reports_host_provenance_mismatch():
    baseline = _report([_config(2)])
    candidate = _report([_config(2)])
    baseline["fixture_sha256"] = candidate["fixture_sha256"] = "a" * 64
    baseline["host"] = {
        "platform": "linux",
        "machine": "x86_64",
        "processor": "cpu-a",
        "python": "3.8.20",
        "cpu_count": 8,
    }
    candidate["host"] = dict(baseline["host"])
    candidate["host"]["cpu_count"] = 16

    result = compare_reports(baseline, candidate)

    assert result["provenance_mismatches"] == ["host.cpu_count"]


def test_compare_release_reports_missing_provenance():
    baseline = _report([_config(2)])
    candidate = _report([_config(2)])
    candidate["fixture_sha256"] = "a" * 64
    candidate["host"] = {
        "platform": "linux",
        "machine": "x86_64",
        "processor": "cpu",
        "python": "3.8.20",
        "cpu_count": 8,
    }

    result = compare_reports(baseline, candidate)

    assert "fixture_sha256" in result["baseline_missing_provenance"]
    assert "host.platform" in result["baseline_missing_provenance"]
    assert result["candidate_missing_provenance"] == []


def test_compare_release_cli_writes_machine_readable_output(tmp_path):
    from benchmarks import compare_release

    baseline = tmp_path / "baseline.json"
    candidate = tmp_path / "candidate.json"
    output = tmp_path / "comparison.json"

    host = {
        "platform": "linux",
        "machine": "x86_64",
        "processor": "cpu",
        "python": "3.8.20",
        "cpu_count": 8,
    }
    baseline.write_text(
        json.dumps({
            "schema_version": 1,
            "fixture_sha256": "a" * 64,
            "host": host,
            "configurations": [_config(2)],
        }),
        encoding="utf-8",
    )
    candidate.write_text(
        json.dumps({
            "schema_version": 1,
            "fixture_sha256": "a" * 64,
            "host": host,
            "configurations": [_config(2)],
        }),
        encoding="utf-8",
    )

    assert compare_release.main([
        str(baseline),
        str(candidate),
        "--require-same-configurations",
        "--require-compatible-environment",
        "--max-regression-percent",
        "10",
        "--output",
        str(output),
    ]) == 0

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["threshold_percent"] == 10.0
    assert payload["require_same_configurations"] is True
    assert payload["require_compatible_environment"] is True
    assert payload["failures"] == []



def test_compare_release_reports_metric_availability_mismatch():
    baseline_config = _config(2)
    candidate_config = _config(2)
    candidate_config.pop("encode_seconds_median")

    result = compare_reports(
        _report([baseline_config]),
        _report([candidate_config]),
    )

    assert result["metric_availability_mismatches"] == [
        (
            (
                2,
                None,
                128,
                "bit_plane",
                2,
                False,
                False,
                False,
            ),
            "encode_seconds_median",
            True,
            False,
        )
    ]


def test_compare_release_zero_higher_is_better_baseline_is_improvement():
    baseline = _report([_config(2, encode_vps=0.0)])
    candidate = _report([_config(2, encode_vps=10.0)])

    result = compare_reports(baseline, candidate, threshold_percent=1.0)
    assert not result["failures"]


def test_compare_release_cli_rejects_metric_availability_mismatch(tmp_path):
    from benchmarks import compare_release

    baseline = tmp_path / "baseline.json"
    candidate = tmp_path / "candidate.json"

    old = _config(2)
    new = _config(2)
    new.pop("decode_seconds_median")

    baseline.write_text(
        json.dumps({"schema_version": 1, "configurations": [old]}),
        encoding="utf-8",
    )
    candidate.write_text(
        json.dumps({"schema_version": 1, "configurations": [new]}),
        encoding="utf-8",
    )

    assert compare_release.main([
        str(baseline),
        str(candidate),
        "--require-complete-metrics",
    ]) == 1
