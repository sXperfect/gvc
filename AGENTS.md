# AGENTS.md — GVC

## Sources of truth

1. `pyproject.toml` — package metadata and dependency policy.
2. `.github/workflows/ci.yml` and `scripts/ci.py` — hosted/local CI.
3. `docs/design/versioning.md` — release-line and Python-floor policy.
4. `docs/design/testing.md` — correctness strategy.
5. `CHANGELOG.md` — user-visible changes.

## GVC 1.0 invariants

- Python support starts at 3.8.
- Patch releases in 1.0.x must not raise that floor.
- A future minor release may deliberately raise the Python floor after the
  previous line remains available as a maintenance/release line.
- Preserve serialized GVC data-structure and codec semantics unless a
  separately documented compatibility change requires otherwise.
- Keep unit tests self-contained: no network, external JBIG executable, or
  large genomic download is required for the core suite.
- Optional VCF parsing and Numba acceleration must not become mandatory imports.
- Do not reintroduce the unmaintained `tspsolve` dependency.
- Avoid removed NumPy aliases and test against newest-resolvable dependencies.

## Verification

Canonical local gate:

```bash
./scripts/verify.sh
```

Useful focused gates:

```bash
python scripts/ci.py metadata
python scripts/ci.py native
python scripts/ci.py test
python scripts/ci.py optional
python scripts/ci.py deps
```

Hosted CI runs only on pushes to main and pull requests targeting main.
