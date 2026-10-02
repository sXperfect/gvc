# CI design

GVC follows the local-first, credit-conscious CI structure used by
`gunz-utils`, adapted to GVC's native Cython/CMake components.

## Hosted triggers

GitHub Actions runs only for:

- pushes to `main`;
- pull requests targeting `main`.

Feature-branch pushes do not trigger hosted CI merely for convenience.
Superseded runs are cancelled and workflow permissions remain read-only.

## Python compatibility coverage

GVC 1.0.x has a Python 3.8 floor. CI uses separate jobs so interpreter failures are isolated:

| Job | Purpose |
| --- | --- |
| Python 3.8 reviewed stack | exact reviewed compatibility ceiling, full automated gates, real JBIG, pinned historical fixture, benchmark smoke |
| Python 3.8 resolver drift | unconstrained eager-resolution probe against the reviewed ceiling |
| Python 3.9-3.14 | compatibility matrix; representative versions run the broader native/packaging gates while all versions exercise maintained tests and dependency checks |

The matrix uses `fail-fast: false` so one interpreter failure does not hide
results from the others. Eager upgrades make dependency ceilings visible rather
than accidentally passing against stale cached packages.

## Local gate

Run:

```bash
./scripts/verify.sh
```

or:

```bash
python scripts/ci.py all
```

The dispatcher checks release/dependency metadata, syntax, Cython and CMake
native builds, pytest, CLI startup, and reports resolved dependency versions.
The `optional` gate verifies that optional integrations import when the
corresponding extras are installed.


## Release validation

Production-scale release validation is intentionally separate from PR CI and is
started only through `workflow_dispatch`. It runs on Linux, uses the reviewed
Python 3.8 stack, verifies the pinned historical fixture, can compare against a
content-addressed controlled-machine benchmark baseline, and emits retained
machine-readable evidence. RC mode additionally builds and inspects the exact
wheel/sdist artifacts intended for promotion.
