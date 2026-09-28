# Changelog

All notable changes to this repository will be documented here.

## Unreleased

### Added

- Explicit GVC 1.0.x Python 3.8 support contract and minor-line Python-floor policy.
- Modern PEP 517/518 package metadata.
- Self-contained pytest coverage for transforms, serialization, solver behavior,
  CLI parsing, and end-to-end codec pipeline round trips.
- Cross-version CI coverage for Python 3.8, 3.9, 3.10, 3.12, and 3.14.
- Dependency-resolution reporting and optional-dependency smoke checks.

### Changed

- GVC 1.0.x uses Cython 3.2.9 as its Python-3.8 build baseline while newer
  interpreters can resolve Cython 3.3+.
- GVC 1.0.x now starts from the newest scientific-stack generation compatible
  with Python 3.8: NumPy 1.24.4 and SciPy 1.10.1.
- Optional VCF, acceleration, and JBIG-example dependencies use modern lower
  bounds while allowing newer Python versions to resolve newer releases.
- The unmaintained `tspsolve` dependency is replaced by the internal
  deterministic nearest-neighbor implementation.
- Deprecated NumPy scalar aliases are removed from maintained code paths.
- VCF and Numba integrations are lazy/optional rather than mandatory imports.
