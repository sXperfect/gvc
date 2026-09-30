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
