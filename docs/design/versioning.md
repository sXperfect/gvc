# Versioning and Python support

GVC uses release lines to decouple modernization from legacy Python support.

## Policy

- The package metadata in `pyproject.toml` is authoritative for the minimum
  Python version of the current release line.
- A minor GVC release may raise the minimum supported Python version.
- Raising the Python floor does not require deleting the older release line.
  Older users can continue installing the newest compatible release.
- Patch releases must not raise the Python floor.
- Major versions are reserved for GVC API, file-format, or codec compatibility
  breaks rather than routine Python lifecycle changes.

## Planned support lines

| GVC line | Minimum Python | Branch policy |
| --- | --- | --- |
| 1.0.x | 3.8 | legacy-compatible baseline |
| 1.1.x | 3.9 | next modernization line |
| 1.2.x | 3.10 | planned future line |

Before `main` moves from one support line to the next, create a maintenance
branch such as `release/1.0`. Backport important fixes there when they remain
applicable.

For example, once GVC 1.1 development starts:

- `release/1.0` retains `requires-python = ">=3.8"`;
- `main` changes to `requires-python = ">=3.9"`;
- CI on each branch tests the Python floor appropriate to that release line.

This lets current development remove obsolete compatibility code without
preventing users on older Python versions from using a maintained historical
release.
