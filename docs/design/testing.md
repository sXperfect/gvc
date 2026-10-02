# Testing strategy

GVC testing is organized around compatibility boundaries rather than historical scripts.

## Layer 1 — pure utilities

Covers bit I/O, permutation serialization, AMax vectors, distance functions, and deterministic nearest-neighbor ordering.

These tests should be fast, deterministic, and free of optional dependencies.

## Layer 2 — transforms

Covers genotype parsing, adaptive missing-value representation, bit-plane binarization, row-binary splitting, inverse transforms, and sorting/unsorting.

Round trips are preferred over implementation-specific assertions.

## Layer 3 — serialized data structures

Parameter sets, blocks, and access units are serialized to bytes and parsed back through the same public structures.

These are compatibility-sensitive regression tests because changes can affect existing .gvc files.

## Layer 4 — codec pipeline

The pipeline is tested end-to-end using a deterministic in-memory lossless matrix codec injected into the historical codec registry.

This deliberately tests GVC-owned logic without requiring an external JBIG executable.

## Layer 5 — integrations

VCF parsing, production JBIG adapters, optional acceleration, and larger datasets belong in explicitly scoped integration tests. They must not make the core unit suite dependent on external tools or downloads.

## Maintained compatibility coverage

The current 1.0.x suite includes:

- byte-exact v1 golden fixtures and structural corpus checks;
- exhaustive truncation rejection and malformed framing/reference tests;
- haploid, diploid, mixed-ploidy, missing-value, phase, and row-bin cases;
- random-access boundary, sample-order, missing-metadata, and invalid-query tests;
- Python/native parity for query expansion, row-bin decoding, permutation
  decoding, and phase reconstruction;
- multiprocessing parity and lifecycle/failure regressions;
- real JBIG and pinned historical compatibility as release gates.

Production-scale datasets, external historical `.gvc` artifacts, and
non-Linux native validation remain controlled/offline release activities.
