# AGENTS.md — GVC

This file defines repository-specific operating rules for automated and human development.

## Sources of truth

1. pyproject.toml — package version, Python floor, dependency groups, test configuration.
2. .github/workflows/ci.yml and scripts/ci.py — hosted/local verification behavior.
3. docs/development/python-support.md — interpreter compatibility policy.
4. CHANGELOG.md — user-visible changes.
5. docs/design/testing.md and docs/development/releases.md — test and release policy.\n6. CONTRIBUTING.md — contributor workflow.\n7. LICENSE and NOTICE.md — licensing and provenance.

Historical material under docs/history/ is preserved for provenance and is not current operating policy.

## Release-line invariants

- GVC 1.x supports Python 3.8 and newer interpreters covered by CI.
- Do not raise the Python minimum within the 1.x line.
- Raising the minimum Python version requires a new major release.
- Python support is encoded in pyproject.toml [project].requires-python.
- The package version in pyproject.toml is the static release source of truth.
- Do not introduce independent hard-coded runtime version literals.

## Compatibility priorities

For the 1.x modernization line, preserve unless an intentional compatibility change is documented and tested:

- serialized GVC data-structure semantics;
- codec identifiers and binarization identifiers;
- decode behavior for existing payloads;
- command-line encode/decode concepts;
- Python 3.8 compatibility.

Prefer additive cleanup over format changes.

## Testing policy

Pytest is the canonical test runner.

Unit tests must be self-contained. Core tests must not require network access, an external JBIG executable, large downloaded genomic datasets, or optional VCF/Numba packages.

For entropy-codec pipeline tests, inject a deterministic in-memory lossless codec into the existing codec contract.

Every bug fix should add a focused regression test where practical.

## CI policy

Hosted CI follows the credit-conscious pattern used by gunz-utils:

- run on pushes to main;
- run on pull requests targeting main;
- do not run merely because a feature branch was pushed;
- keep one primary verify job;
- use minimal contents: read permissions;
- cancel superseded runs;
- mirror hosted checks locally through scripts/ci.py.

Canonical local gate:

    bash scripts/verify.sh

## Dependency policy

Keep the core installation as small as practical.

- NumPy and SciPy are core numerical dependencies.
- cyvcf2 belongs to the vcf extra.
- Numba belongs to the speed extra; maintained code must have a correct non-Numba path.
- Do not reintroduce tspsolve; nearest-neighbor ordering is maintained internally.
- Build-only dependencies belong in [build-system].

Before adding a dependency, prefer a small local implementation when the required behavior is narrow and can be comprehensively tested.

## Native code

Cython extensions and library/libgvc are compatibility-sensitive.

- Avoid -march=native in portable builds.
- Keep compiler optimizations architecture-neutral.
- Run the native build gate after changes to Cython, C, CMake, or packaging.
- Native acceleration must preserve observable Python semantics.

## Definition of done

Applicable work is complete only when Python 3.8 support is preserved for GVC 1.x, relevant unit/regression tests exist, the consolidated gate passes in a prepared environment, dependency boundaries remain intentional, documentation matches behavior, provenance is preserved, and user-visible changes are reflected in CHANGELOG.md.
