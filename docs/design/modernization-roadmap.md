# GVC modernization roadmap

## Phase 1 — provenance and compatibility foundation

- import the historical implementation with its BSD license;
- establish GVC 1.x as Python >=3.8;
- document that a higher Python floor requires a major release;
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

## Phase 5 — stabilize GVC 1.x

- complete documentation and installation guidance;
- establish a release branch when the first maintained 1.x release is cut;
- keep Python 3.8 compatibility fixes isolated from future major-line work.

## Phase 6 — future GVC 2.x

Only after the 1.x baseline is reproducible and audited:

- choose a newer Python minimum;
- modernize dependencies around that floor;
- consider stronger typing/lint gates;
- consider format/API changes that require a major version.
