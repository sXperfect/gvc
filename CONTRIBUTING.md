# Contributing to GVC

GVC 1.x is the compatibility and modernization line. Changes should improve correctness, reproducibility, maintainability, or documentation without silently raising the Python floor or changing the serialized codec contract.

## Development setup

Python 3.8 is the minimum supported interpreter for the 1.x line.

Install the core project and test dependencies:

    python -m pip install -e ".[test]"

Optional VCF integration:

    python -m pip install -e ".[test,vcf]"

Optional Numba acceleration:

    python -m pip install -e ".[test,speed]"

## Verification

Run the canonical local gate:

    bash scripts/verify.sh

Focused gates:

    python scripts/ci.py metadata
    python scripts/ci.py syntax
    python scripts/ci.py native
    python scripts/ci.py test

## Test expectations

Prefer small deterministic unit/regression fixtures.

Core tests must not depend on network access, downloaded genomic datasets, or an external JBIG executable. Test entropy-codec integration by injecting a lossless in-memory codec through GVC's codec contract.

For changes to serialized structures, add binary round-trip tests. For numerical or algorithmic changes, test edge cases and invariants before optimizing the implementation.

## Compatibility

Within GVC 1.x:

- keep Python >=3.8;
- preserve codec/binarization identifiers;
- preserve existing binary semantics unless a documented compatibility fix requires otherwise;
- do not introduce a dependency that unnecessarily narrows supported Python versions.

A deliberate minimum-Python increase belongs in the next major release.

## Pull requests

Pull requests should state:

- the behavior being changed;
- compatibility implications;
- verification performed;
- whether serialized data structures or codec behavior are affected.

Keep unrelated modernization work separate when practical.
