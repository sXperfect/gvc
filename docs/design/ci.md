# CI design

GVC adopts the useful structural principles of the gunz-utils CI without copying gates that do not yet fit this codebase.

## Goals

- local and hosted verification use the same commands;
- CI failures identify real regressions rather than known historical setup gaps;
- Python 3.8 remains continuously verified for the 1.x line;
- native Cython/C code is built in CI;
- unit tests do not depend on external JBIG binaries or genomic downloads;
- feature-branch pushes do not consume hosted CI unnecessarily.

## Hosted workflow

The single primary job is CI / verify.

It runs for pushes to main and pull requests targeting main. The job uses minimal read permissions, a timeout, and concurrency cancellation.

The initial compatibility sequence is:

1. Python 3.8: install core + test dependencies;
2. syntax/bytecode compilation;
3. Cython import/build verification;
4. standalone CMake helper build;
5. full pytest suite;
6. Python 3.10: reinstall and run the full pytest suite as a newer-interpreter compatibility check.

## Local dispatcher

scripts/ci.py is the canonical dispatcher. scripts/verify.sh invokes its all gate.

The first modernization phase intentionally does not make Ruff or mypy blocking gates. They can be introduced after the imported legacy code is cleaned enough that those checks represent regressions rather than an inventory of pre-existing style/type debt.
