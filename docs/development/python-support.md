# Python and dependency support for GVC 1.0.x

GVC 1.0.x intentionally keeps a broad compatibility window while allowing
dependencies to advance as far as each Python interpreter permits.

## Supported interpreters

The tested 1.0.x window is Python **3.8 through 3.12**. Package metadata uses:

```text
requires-python = ">=3.8,<3.13"
```

Python 3.13+ is intentionally left to a later GVC line until the native/Cython
surface and NumPy 2.x behavior are validated explicitly.

## Dependency strategy

Dependencies use environment markers instead of one historical pin.

- Python 3.8 remains on the newest compatible NumPy 1.24 / SciPy 1.10 family.
- Python 3.9-3.12 may resolve newer NumPy 1.x and SciPy 1.x releases.
- `cyvcf2` is optional because only VCF I/O requires it.
- Numba is optional because GVC has a correct non-Numba path.
- `tspsolve` has been replaced by the maintained in-repository deterministic
  nearest-neighbor implementation.
- Pillow is not a core dependency; it is only referenced by historical external
  JBIG integration documentation.

CI installs the `all,test` extras on every supported Python version so both
the core and optional dependency ranges are continuously exercised.

## Patch-release rule

Within 1.0.x, dependency range changes are acceptable when they preserve the
documented Python window, serialized format behavior, and public API behavior.
A dependency update that requires dropping a supported Python version belongs
in the next minor release line instead.
