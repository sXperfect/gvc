# GVC benchmarks

This directory contains reproducible release-hardening benchmarks. They are
intended to establish and compare baselines; ordinary CI does **not** enforce
wall-clock thresholds because GitHub-hosted runner timing is too noisy for a
stable performance contract.

Each worker-count configuration is executed in a fresh Python subprocess so
peak RSS and child-process RSS are isolated rather than accumulated across the
whole benchmark session.

## Metrics

The JSON report records:

- encode and decode wall time;
- encode/decode variants per second;
- full-sample random-access latency for one sample;
- GVC stream bytes and metadata-sidecar bytes;
- decoded genotype-text bytes;
- decoded-genotype-bytes / total-GVC-bytes size ratio;
- peak RSS for the benchmark process and its children;
- Python/platform information;
- raw per-repetition samples plus median/min/max summaries.

The size ratio is a genotype-text storage comparison, **not** a VCF compression
ratio. The input VCF contains annotations that GVC does not encode, so comparing
the whole VCF byte size directly with a genotype-only GVC stream would be
misleading.

## Example

With JBIG-KIT installed:

```bash
python benchmarks/run_release.py \
  tmp/historical/test_block01.vcf.gz \
  --workers 0 1 2 4 8 \
  --repetitions 5 \
  --block-size 2048 \
  --output benchmarks/results/linux-py38.json
```

For row-bin-split:

```bash
python benchmarks/run_release.py \
  tmp/historical/test_block01.vcf.gz \
  --workers 0 1 2 4 \
  --repetitions 3 \
  --block-size 2048 \
  --binarization row_bin_split \
  --axis 0 \
  --output benchmarks/results/row-bin-split.json
```

Use `--start-method spawn` when validating Windows/macOS-like process
semantics on Linux.

## Baseline policy

Do not commit arbitrary hosted-runner timing thresholds as pass/fail tests.
Before a release candidate, capture at least one controlled-machine baseline
and retain the JSON report externally or under a release-specific benchmark
record. Compare later runs on the same machine/toolchain and investigate
material regressions in throughput, memory, or storage ratio.

Production release validation should additionally use the full LUH historical
fixture, realistic worker counts, sorting enabled where relevant, and the real
JBIG-KIT executables.


## Comparing controlled-machine baselines

After capturing two reports on the same machine/toolchain, compare them with:

```bash
python benchmarks/compare_release.py \
  benchmarks/results/v1.0.1rc1-baseline.json \
  benchmarks/results/candidate.json
```

This is report-only by default. For a controlled release gate, an explicit
regression budget can be supplied:

```bash
python benchmarks/compare_release.py \
  baseline.json candidate.json \
  --require-same-configurations \
  --max-regression-percent 10
```

Throughput metrics treat lower values as regressions; latency, RSS, and encoded
size treat higher values as regressions. Do not use these thresholds on hosted
CI runners unless the environment is demonstrably stable.


## Canonical controlled baseline

A retained release baseline is only considered directly comparable when all of
the following hold:

- the fixture SHA-256 matches exactly;
- benchmark configuration keys match exactly;
- platform, machine architecture, processor identity, Python version, and CPU
  count match;
- the run was captured on a controlled machine rather than a shared hosted
  runner;
- the baseline JSON itself is retained immutably and identified by SHA-256.

For RC/final validation, supplying a baseline also requires an explicit
`--max-regression-percent` budget. This prevents a baseline comparison from
being treated as a release gate without an agreed threshold.

The candidate and retained baseline may come from different GVC commits or
package versions—that is the purpose of the comparison—but differences in
fixture or controlled-machine provenance must be treated as non-comparable
rather than as performance regressions.
