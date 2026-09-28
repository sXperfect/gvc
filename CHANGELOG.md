# Changelog

All notable changes to this repository will be documented here.

## Unreleased

### Added

- Explicit GVC 1.x Python support policy with Python 3.8 as the minimum.
- Modern PEP 517/518 package metadata.
- Self-contained pytest-based verification.
- Local-first CI dispatcher and GitHub Actions workflow.
- Provenance and modernization design documentation.

### Changed

- VCF parsing and Numba acceleration are optional dependency surfaces.
- Legacy dependency pins are replaced by compatibility ranges suitable for the Python 3.8 release line.
- Nearest-neighbor ordering no longer requires the unmaintained tspsolve package.
- Deprecated NumPy scalar aliases are removed from maintained paths.

### Notes

No GVC 2.x compatibility-floor change is part of this work. A higher minimum Python version will be considered only in a future major release.
