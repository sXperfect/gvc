# Python and dependency support for GVC 1.0.x

GVC 1.0.x has a minimum supported interpreter of **Python 3.8**.

The package intentionally has no artificial upper Python bound. CI currently
tests CPython 3.8 through 3.14. When a newer Python release appears, it is
treated as a compatibility target for the existing 1.0.x line unless a real
implementation or dependency incompatibility is demonstrated.

## Release-line rule

- 1.0.x keeps Python >=3.8.
- Patch releases must not raise that floor.
- A future minor line may deliberately raise the floor, for example 1.1.x to
  Python >=3.9, while 1.0.x remains installable for Python 3.8 users.
- Major releases are reserved for API, file-format, or codec compatibility
  changes rather than routine interpreter lifecycle updates.

## Dependency strategy

GVC declares the newest dependency generation that still supports Python 3.8
as its lower bound, plus a next-major upper bound:

- NumPy >=1.24.4,<3
- SciPy >=1.10.1,<2
- cyvcf2 0.31.4 on Python 3.8 (last release with CPython 3.8 wheels); cyvcf2 >=0.34.0,<1 on Python >=3.9
- Numba >=0.58.1,<1 (optional acceleration)
- Pillow >=10.4.0,<13 (optional JBIG integration support)
- pytest >=8.3.5,<10
- Cython >=3.2.9,<4
- reviewed Python 3.8 build tooling: setuptools 75.3.4, wheel 0.45.1, and build 1.2.2

pip then uses each project's Requires-Python metadata to select the newest
compatible release on the active interpreter. Thus Python 3.8 retains a modern
but compatible dependency generation, while newer interpreters exercise newer
NumPy, SciPy, Cython, pytest, Numba, cyvcf2, and Pillow releases. The VCF extra uses an explicit interpreter marker because cyvcf2 0.33.0 declares Python 3.8 compatibility but no longer publishes CPython 3.8 wheels.

CI uses eager upgrades deliberately to make these dependency ceilings visible
rather than accidentally passing against stale cached packages.


## Python 3.8 compatibility gate

Python 3.8 is treated as the compatibility anchor for the 1.0.x line.

`ci/constraints/py38-latest.txt` records the reviewed newest stack currently
known to support Python 3.8. CI uses it in two complementary ways:

1. **reviewed stack** — install the exact versions and run the complete native,
   serialization, golden-format, pipeline, VCF, CLI, and optional-integration
   suite;
2. **resolver probe** — create a second clean Python 3.8 environment without
   constraints, request eager upgrades, and compare what pip resolves against
   the reviewed ceiling.

If an upstream project publishes a newer Python-3.8-compatible release, the
resolver probe fails intentionally. The new version should be tested before the
constraint is updated. This prevents both accidental stagnation and unreviewed
dependency drift in the legacy-compatible release line.
