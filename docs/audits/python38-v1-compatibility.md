# GVC 1.0.x / Python 3.8 compatibility audit

## Scope

This audit treats CPython 3.8 as the compatibility anchor for the GVC 1.0.x
maintenance line. The goal is to modernize dependencies and implementation
without silently changing the historical GVC serialization contract.

## Implemented verification

The maintained test suite now covers:

- reviewed latest Python 3.8 dependency stack plus an unconstrained resolver
  probe;
- byte-exact v1 golden serialization for parameter sets, permutations, AMax,
  blocks, access units, and a complete structural .gvc fixture;
- exhaustive truncation rejection for the structural golden fixture;
- malformed sizes, IDs, flags, payload lengths, permutation values, and
  bitstream boundaries;
- deterministic native-vs-Python differential tests for cquery, cdebinarize,
  crc_id, and the standalone libgvc permutation helper;
- generated encode/decode cases across ploidy 1-4, phase modes, missing/NA
  values, sorting, transpose, and both binarization schemes;
- tiny diploid and haploid VCF fixtures, including mixed phasing,
  multiallelic calls, missing genotypes, partial final blocks, and block
  boundaries;
- complete Encoder -> .gvc -> Decoder file round trips for bit-plane,
  row-bin-split, and haploid data using a deterministic JBIG-header-compatible
  in-memory codec;
- wheel and sdist build/install isolation from outside the source checkout;
- installed native-extension imports and installed CLI startup.

## Correctness fixes found by the audit

The expanded tests exposed and now protect fixes for:

- cyvcf2 phasing convention conversion at the ingestion boundary;
- haploid zero-width phase matrices;
- row-bin-split tensor row-count reconstruction;
- Encoder's default codec name not matching the registered codec;
- unsafe standalone C permutation decoding, including floor(log2(n)) for
  non-power-of-two permutation sizes and unchecked payload reads;
- unknown sample IDs silently mapping to sample column zero;
- interval queries before the first block or inside metadata gaps;
- BinMat emitting dimensions but omitting the matrix payload;
- truncated lazy payload regions being accepted through seek arithmetic;
- explicit serialized integer/flag/payload size validation.

## Hosted CI policy

Python 3.8 runs the complete release-quality gate, including packaging
isolation. Newer interpreters exercise compatibility, native modules, tests,
CLI, and optional integrations without rebuilding release artifacts repeatedly.
This keeps CI credit-conscious while preserving the strongest checks at the
supported floor.

## Offline / external verification still required

The following checks cannot be fully represented by the self-contained hosted
suite and should be performed before declaring a production 1.0.x release:

1. Decode representative historical .gvc files produced by the original LUH
   implementation, not only synthetic structural golden fixtures.
2. Run encode/decode with the actual chosen external JBIG executable and verify
   byte/shape interoperability. The hosted tests intentionally inject an
   in-memory codec and do not validate a third-party JBIG binary.
3. Run the historical large VCF fixtures (including test_block01.vcf.gz) and
   compare reconstructed genotypes against the source data.
4. Exercise random-access queries against real metadata sidecars for genomic
   intervals and sample subsets at realistic scale.
5. Exercise the multiprocessing encoder path with multiple workers and verify
   deterministic block ordering and clean worker termination.
6. Stress very large matrices/block counts near serialized field-size
   boundaries and monitor memory use.
7. Run native builds on any production platforms beyond the Linux CI target,
   especially compiler/OpenMP and shared-library behavior.
8. Record compression ratio and throughput regressions against a known 1.0
   baseline once the correctness line is frozen.

These are release checks, not reasons to weaken or skip the self-contained CI
gate.
