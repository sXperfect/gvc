# CI design

The CI design follows the local-first structure used by `gunz-utils`, adapted
to GVC's native extensions and legacy dependency surface.

## Hosted policy

GitHub Actions runs only for pushes to `main` and pull requests targeting
`main`. Feature-branch pushes do not consume hosted CI runs.

The hosted workflow uses one `verify` job with minimal read-only permissions,
a timeout, and concurrency cancellation. GVC 1.0 validates Python 3.8 first and
then Python 3.9 in the same job.

## Local gate

Run:

```bash
./scripts/verify.sh
```

or:

```bash
python scripts/ci.py all
```

The current gates are:

1. package version and Python-floor metadata;
2. Python syntax compilation;
3. Cython extension and CMake native-library builds;
4. the unittest suite;
5. CLI import/help smoke testing.

The older upstream repository stored very large VCF fixtures. The migrated
repository does not require those datasets for its baseline CI: file-backed
integration tests skip when the legacy fixture is absent, while in-memory
round-trip tests install a lossless test-only matrix codec so core encode/decode
behavior remains exercised without an external JBIG executable.

Linting, typing, packaging isolation, and documentation gates should be added
after the legacy implementation is brought to a clean baseline rather than
making the initial CI permanently red for unrelated historical debt.
