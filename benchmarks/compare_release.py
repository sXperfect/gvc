#!/usr/bin/env python3
"""Compare two GVC release-benchmark JSON reports."""

import argparse
import json
import math
import sys
from pathlib import Path


KEY_FIELDS = (
    "workers",
    "start_method",
    "block_size",
    "binarization",
    "axis",
    "sort_rows",
    "sort_cols",
    "transpose",
)

METRICS = {
    "encode_variants_per_second_median": "higher",
    "decode_variants_per_second_median": "higher",
    "encode_seconds_median": "lower",
    "decode_seconds_median": "lower",
    "random_access_seconds_median": "lower",
    "peak_rss_kib_self_median": "lower",
    "peak_rss_kib_children_median": "lower",
    "total_gvc_bytes": "lower",
}


def _load(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError(
            "unsupported benchmark schema_version: {!r}".format(
                data.get("schema_version")
            )
        )
    configurations = data.get("configurations")
    if not isinstance(configurations, list):
        raise ValueError("benchmark report must contain a configurations list")
    return data


def _key(config):
    return tuple(config.get(field) for field in KEY_FIELDS)


def _index(report):
    indexed = {}
    for config in report["configurations"]:
        key = _key(config)
        if key in indexed:
            raise ValueError("duplicate benchmark configuration: {!r}".format(key))
        indexed[key] = config
    return indexed


def _change_percent(baseline, candidate, direction):
    baseline = float(baseline)
    candidate = float(candidate)
    if not math.isfinite(baseline) or not math.isfinite(candidate):
        raise ValueError("benchmark metrics must be finite")
    if baseline == 0:
        if candidate == 0:
            return 0.0
        if direction == "higher":
            return -math.inf
        return math.inf

    if direction == "higher":
        return (baseline - candidate) / baseline * 100.0
    return (candidate - baseline) / baseline * 100.0


PROVENANCE_FIELDS = ("fixture_sha256",)
HOST_PROVENANCE_FIELDS = ("platform", "machine", "processor", "python", "cpu_count")


def _provenance_mismatches(baseline, candidate):
    mismatches = []
    for field in PROVENANCE_FIELDS:
        old = baseline.get(field)
        new = candidate.get(field)
        if old is not None and new is not None and old != new:
            mismatches.append(field)

    baseline_host = baseline.get("host") or {}
    candidate_host = candidate.get("host") or {}
    for field in HOST_PROVENANCE_FIELDS:
        old = baseline_host.get(field)
        new = candidate_host.get(field)
        if old is not None and new is not None and old != new:
            mismatches.append("host." + field)
    return mismatches


def _missing_provenance(report):
    missing = []
    for field in PROVENANCE_FIELDS:
        if report.get(field) is None:
            missing.append(field)
    host = report.get("host") or {}
    for field in HOST_PROVENANCE_FIELDS:
        if host.get(field) is None:
            missing.append("host." + field)
    return missing


def compare_reports(baseline, candidate, threshold_percent=None):
    baseline_index = _index(baseline)
    candidate_index = _index(candidate)

    missing = sorted(set(baseline_index) - set(candidate_index), key=repr)
    extra = sorted(set(candidate_index) - set(baseline_index), key=repr)

    rows = []
    failures = []
    metric_availability_mismatches = []
    for key in sorted(set(baseline_index) & set(candidate_index), key=repr):
        old = baseline_index[key]
        new = candidate_index[key]
        metric_rows = []
        for metric, direction in METRICS.items():
            old_present = metric in old and old[metric] is not None
            new_present = metric in new and new[metric] is not None
            if old_present != new_present:
                metric_availability_mismatches.append(
                    (key, metric, old_present, new_present)
                )
                continue
            if not old_present:
                continue

            regression = _change_percent(old[metric], new[metric], direction)
            record = {
                "metric": metric,
                "baseline": old[metric],
                "candidate": new[metric],
                "regression_percent": regression,
            }
            metric_rows.append(record)
            if (
                threshold_percent is not None
                and regression > threshold_percent
            ):
                failures.append((key, record))

        rows.append({"configuration": key, "metrics": metric_rows})

    return {
        "missing_configurations": missing,
        "extra_configurations": extra,
        "provenance_mismatches": _provenance_mismatches(baseline, candidate),
        "baseline_missing_provenance": _missing_provenance(baseline),
        "candidate_missing_provenance": _missing_provenance(candidate),
        "comparisons": rows,
        "metric_availability_mismatches": metric_availability_mismatches,
        "failures": failures,
    }


def _format_key(key):
    return ", ".join(
        "{}={}".format(field, value)
        for field, value in zip(KEY_FIELDS, key)
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compare equivalent GVC release benchmark configurations."
    )
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument(
        "--max-regression-percent",
        type=float,
        help=(
            "fail when any comparable metric regresses by more than this "
            "percentage; omit for report-only mode"
        ),
    )
    parser.add_argument(
        "--require-same-configurations",
        action="store_true",
        help="fail if either report has unmatched configurations",
    )
    parser.add_argument(
        "--require-compatible-environment",
        action="store_true",
        help="fail when fixture or benchmark host provenance differs",
    )
    parser.add_argument(
        "--require-complete-metrics",
        action="store_true",
        help="fail if comparable metric availability differs between reports",
    )
    parser.add_argument(
        "--output",
        help="optional JSON file for machine-readable comparison results",
    )
    args = parser.parse_args(argv)

    if (
        args.max_regression_percent is not None
        and args.max_regression_percent < 0
    ):
        parser.error("--max-regression-percent must be non-negative")

    try:
        result = compare_reports(
            _load(args.baseline),
            _load(args.candidate),
            args.max_regression_percent,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    result["threshold_percent"] = args.max_regression_percent
    result["require_same_configurations"] = bool(args.require_same_configurations)
    result["require_compatible_environment"] = bool(args.require_compatible_environment)
    result["require_complete_metrics"] = bool(args.require_complete_metrics)

    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True, default=list) + "\n",
            encoding="utf-8",
        )

    for row in result["comparisons"]:
        print(_format_key(row["configuration"]))
        for metric in row["metrics"]:
            print(
                "  {:40s} {:>12} -> {:>12}  regression={:+.2f}%".format(
                    metric["metric"],
                    str(metric["baseline"]),
                    str(metric["candidate"]),
                    metric["regression_percent"],
                )
            )

    if result["provenance_mismatches"]:
        print("benchmark provenance differs:", file=sys.stderr)
        for field in result["provenance_mismatches"]:
            print("  " + field, file=sys.stderr)
    if result["baseline_missing_provenance"]:
        print("baseline missing provenance:", file=sys.stderr)
        for field in result["baseline_missing_provenance"]:
            print("  " + field, file=sys.stderr)
    if result["candidate_missing_provenance"]:
        print("candidate missing provenance:", file=sys.stderr)
        for field in result["candidate_missing_provenance"]:
            print("  " + field, file=sys.stderr)

    if result["metric_availability_mismatches"]:
        print("benchmark metric availability differs:", file=sys.stderr)
        for key, metric, old_present, new_present in result["metric_availability_mismatches"]:
            print(
                "  {}: {} baseline_present={} candidate_present={}".format(
                    _format_key(key),
                    metric,
                    old_present,
                    new_present,
                ),
                file=sys.stderr,
            )

    if result["missing_configurations"]:
        print("missing candidate configurations:", file=sys.stderr)
        for key in result["missing_configurations"]:
            print("  " + _format_key(key), file=sys.stderr)
    if result["extra_configurations"]:
        print("extra candidate configurations:", file=sys.stderr)
        for key in result["extra_configurations"]:
            print("  " + _format_key(key), file=sys.stderr)

    if result["failures"]:
        print("benchmark regressions exceeded threshold:", file=sys.stderr)
        for key, metric in result["failures"]:
            print(
                "  {}: {} {:+.2f}%".format(
                    _format_key(key),
                    metric["metric"],
                    metric["regression_percent"],
                ),
                file=sys.stderr,
            )
        return 1

    if args.require_same_configurations and (
        result["missing_configurations"] or result["extra_configurations"]
    ):
        return 1
    if args.require_compatible_environment and (
        result["provenance_mismatches"]
        or result["baseline_missing_provenance"]
        or result["candidate_missing_provenance"]
    ):
        return 1
    if args.require_complete_metrics and result["metric_availability_mismatches"]:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
