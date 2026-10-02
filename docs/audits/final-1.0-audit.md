# Final GVC 1.0.x pre-release audit

Audit branch: `feat/ci-versioning`

This audit records the repository state immediately before creation of the
`release/1.0` maintenance branch.

## Automated release foundation

The 1.0.x line now has automated coverage for:

- Python 3.8 reviewed compatibility stack;
- Python 3.8 unconstrained resolver drift;
- Python 3.9 through 3.14 compatibility;
- native Cython and standalone CMake components;
- wheel/sdist construction, structural inspection, isolated installation, and
  installed CLI/native imports;
- optional dependency imports;
- byte-exact v1 golden serialization and malformed-stream rejection;
- historical LUH fixture compatibility;
- real JBIG-KIT integration;
- multiprocessing parity, lifecycle, cancellation, abrupt-exit, SIGTERM,
  partial-start, and transactional-output behavior;
- random-access correctness;
- benchmark smoke generation and controlled baseline comparison tooling.

## Release-policy invariants

The maintained 1.0.x branch must preserve:

- Python `>=3.8`;
- the v1 serialized format;
- Linux as the validated release platform;
- macOS/Windows as best-effort unless separate validation evidence is added;
- manual-only production release validation;
- exact package/tag agreement for RC/final versions.

## Packaging and provenance

Release tooling records or verifies:

- package version from `gvc/_version.py`;
- wheel metadata and native extension presence;
- source-distribution build inputs;
- artifact SHA-256 hashes;
- Git commit and clean-tree state;
- historical fixture SHA-256 and Git blob identity;
- benchmark fixture/environment provenance;
- benchmark comparison evidence and explicit regression budget.

## Branch state

Pull request #2 (`feat/ci-versioning` -> `main`) is reported by GitHub as
mergeable with a clean merge state at the time of this audit.

## Remaining release-process gates

These are deliberate next-stage actions rather than unresolved implementation
defects:

1. obtain a final green CI run for the current audit/preflight head;
2. merge/freeze the approved preparation state and create `release/1.0`;
3. promote `1.0.1.dev0` to `1.0.1rc1`;
4. update the changelog for the RC;
5. run the manual controlled release-validation workflow in RC mode using the
   full historical fixture and sorting coverage;
6. retain and review the generated evidence, artifacts, benchmark comparison,
   and hashes;
7. tag the exact reviewed RC commit only after evidence acceptance.

## Audit conclusion

No known code, packaging, compatibility, or CI architecture issue remains that
must be redesigned before creating the 1.0 maintenance branch. The remaining
work is release execution and controlled/offline validation.
