# Versioning and Python support

GVC uses release lines to decouple modernization from legacy Python support.

## Python-floor policy

- GVC 1.0.x requires Python 3.8 or newer.
- Patch releases within 1.0.x do **not** raise the Python floor.
- A new minor line may raise the floor when that change is explicit and
  documented. For example, GVC 1.1.x may require Python 3.9 while 1.0.x
  remains available to Python 3.8 users.
- Major versions remain reserved for GVC API, serialized-format, or codec
  compatibility breaks rather than routine Python lifecycle updates.

Planned progression:

| GVC line | Minimum Python |
| --- | ---: |
| 1.0.x | 3.8 |
| 1.1.x | 3.9 |
| 1.2.x | 3.10 |

Before main moves to a new minor line, create a maintenance branch such as
`release/1.0` when continued fixes for the previous Python floor are desired.

## Dependency policy for 1.0.x

The 1.0.x package uses the newest dependency generation that can still be
installed on Python 3.8 as its lower bound:

- NumPy >= 1.24.4
- SciPy >= 1.10.1
- cyvcf2 >= 0.33.0 for the optional VCF surface
- Numba >= 0.58.1 for optional acceleration
- Pillow >= 10.4.0 for the optional JBIG integration example
- pytest >= 8.3.5 for tests
- Cython >= 3.2.9 for builds/tests

There are intentionally no global upper pins for these runtime/test
dependencies. pip uses each dependency's own `Requires-Python` metadata to
choose the newest resolvable release for the active interpreter. This means a
Python 3.8 environment receives the newest compatible legacy generation while
newer Python environments exercise newer NumPy/SciPy and toolchain releases.

CI therefore tests both dimensions: the Python 3.8 compatibility floor and
representative newer interpreters with eager dependency upgrades.
