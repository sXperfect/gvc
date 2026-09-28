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

GVC 1.0.x has a Python 3.8 floor. CI sequentially exercises representative
interpreter/dependency transitions in one hosted job:

| Python | Purpose |
| --- | --- |
| 3.8 | minimum interpreter + newest dependencies still resolvable there |
| 3.9 | first post-floor dependency transition |
| 3.10 | full native/core/optional gate |
| 3.11 | intermediate compatibility regression coverage |
| 3.12 | full current scientific-stack gate |
| 3.13 | intermediate compatibility regression coverage |
| 3.14 | forward-compatible core/native gate |

Each environment uses pip's eager upgrade strategy so CI does not accidentally
pass against stale cached dependencies.

Python 3.8 first runs the core installation without optional integrations and
then installs all optional extras. Python 3.9-3.12 run the full optional stack;
Python 3.14 runs the core surface so a lagging optional package does not
artificially redefine GVC's core interpreter compatibility.

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
