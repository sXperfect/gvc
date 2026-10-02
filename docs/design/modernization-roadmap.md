# GVC modernization roadmap

## Phase 1 — provenance and compatibility foundation

- import the historical implementation with its BSD license;
- establish GVC 1.x as Python >=3.8;
- document release-line Python floors so a future minor line may raise the
  minimum while 1.0.x remains available for Python 3.8;
- replace legacy packaging metadata with pyproject.toml.

## Phase 2 — test foundation

- replace malformed/external-fixture legacy tests with pytest;
- test bitstreams, binarization, sorting, payloads, and decode round trips;
- isolate third-party entropy codecs behind test doubles;
- keep unit tests deterministic and network-free.

## Phase 3 — dependency and native-build modernization

- remove obsolete tspsolve by maintaining the small nearest-neighbor routine locally;
- make VCF parsing and Numba acceleration optional;
- remove deprecated NumPy aliases;
- make Cython compilation portable and keep CMake builds warning-aware.

## Phase 4 — correctness audit

- add regression tests for each serialized data structure;
- audit boundary cases for empty/singleton matrices, missing genotypes, ploidy, row/column permutations, and random access;
- compare Python and accelerated implementations;
- add small benchmark/regression fixtures only after correctness gates are stable.

## Phase 5 — stabilize GVC 1.0.x

- complete release documentation and installation guidance;
- establish `release/1.0` when the maintained 1.0.x line is frozen;
- retain Python >=3.8 and the v1 serialized format for all 1.0.x patches.

## Phase 6 — future minor/major lines

After the 1.0.x baseline is reproducible and audited:

- a future minor line may deliberately raise the Python minimum;
- modernize dependencies around that new floor without rewriting 1.0.x history;
- reserve major versions for breaking API, serialized-format, or codec changes;
- consider stronger typing/lint gates independently of format compatibility.
