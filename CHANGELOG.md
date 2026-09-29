# Changelog

All notable changes to this repository will be documented here.

## Unreleased

### Added

- Native Cython extension regression tests for query-index expansion,
  row-split decoding, and permutation decoding across the Python/NumPy matrix.
- Byte-exact v1 golden-format tests, exhaustive byte-truncation rejection,
  malformed-stream boundaries, tiny VCF integration fixtures, and generated
  ploidy/phase pipeline cases.

- Explicit GVC 1.0.x Python 3.8 support contract and minor-line Python-floor policy.
- Modern PEP 517/518 package metadata.
- Native CI verification now exercises PEP-517-built Cython modules instead of invoking `setup.py` directly.
- Self-contained pytest coverage for transforms, serialization, solver behavior,
  CLI parsing, and end-to-end codec pipeline round trips.
- Cross-version CI coverage for Python 3.8 through 3.14, with Python 3.8 as the reviewed compatibility anchor.
- Dependency-resolution reporting and optional-dependency smoke checks.
- Python 3.8 wheel and sdist build/install smoke tests from outside the source checkout.
- Python 3.8 wheel/sdist build-and-install isolation gate with native-extension and installed-CLI smoke tests.
- Explicit source-distribution manifest for the Cython `.pyx` build sources.
- Additional byte-exact v1 structural fixtures for row-bin-split/missing values,
  mixed phase payloads, and sorted multi-plane payloads.
- End-to-end random-access and multiprocessing parity regressions.
- A parent-supervised multiprocessing pipeline with structured errors,
  bounded backpressure, progress/watchdog support, transactional output,
  failure cleanup, and fork/spawn lifecycle coverage.
- A maintained JBIG-KIT subprocess integration with executable discovery,
  timeout/error validation, and real codec release-gate tests.
- A pinned historical LUH VCF compatibility gate sourced from upstream commit
  f9af2127a2ff0b87727924860e33f1905fd507cd.
- An isolated release benchmark harness for encode/decode throughput, peak RSS,
  storage size, and random-access latency.
- Bounded real-JBIG multiprocessing queue-pressure and random-access release
  regressions.
- Complete file-level Encoder/Decoder regression coverage for bit-plane,
  row-bin-split, and haploid VCF inputs.
- Random-access index, BinMat framing, and standalone libgvc safety tests.

### Changed

- The Python 3.8 resolver probe promoted setuptools 75.3.4 as the reviewed
  build-tool ceiling after detecting it as the newest compatible release.
- GVC 1.0.x uses Cython 3.2.9 as its Python-3.8 build baseline while newer
  interpreters can resolve Cython 3.3+.
- GVC 1.0.x now starts from the newest scientific-stack generation compatible
  with Python 3.8: NumPy 1.24.4 and SciPy 1.10.1.
- Optional VCF, acceleration, and JBIG-example dependencies use modern lower
  bounds while allowing newer Python versions to resolve newer releases.
- Python 3.8 VCF integration is pinned to the last wheel-backed cyvcf2 0.31.x
  line; Python 3.9+ uses cyvcf2 0.34+.
- The unmaintained `tspsolve` dependency is replaced by the internal
  deterministic nearest-neighbor implementation.
- Deprecated NumPy scalar aliases are removed from maintained code paths.
- Phasing text reconstruction now consistently follows the GVC convention
  `0 = |`, `1 = /` in Python and native helpers.
- Native Cython helpers validate malformed dimensions before entering
  bounds-check-disabled loops.
- Lazy genotype payload regions now validate physical byte availability before
  seeking, so truncated files cannot be accepted by length arithmetic alone.
- Haploid genotype text reconstruction now uses an indexable one-allele
  codebook.
- Unknown sample IDs now fail explicitly instead of silently mapping to column zero.
- Row-bin-split shape reconstruction and haploid phase handling are corrected.
- `BinMat` now writes and reads its matrix payload rather than dimensions only.
- The standalone C permutation bridge uses a length-aware checked decoder.
- VCF ingestion now finalizes metadata on exact block boundaries and splits
  blocks when ploidy changes, preserving one ploidy per ParameterSet.
- Random access handles sample-only queries and empty intervals deterministically.
- Spawned encoder workers can rebuild process-local codec/plugin state through
  an explicit picklable initializer hook.
- The standalone ctypes/C permutation decoder now rejects duplicate IDs and
  trailing payload bytes just like the Python reference implementation.
- VCF and Numba integrations are lazy/optional rather than mandatory imports.
