# Release process

## Version model

GVC uses semantic release lines together with an explicit Python support policy.

The maintained 1.x line has:

- package version source: pyproject.toml;
- minimum Python: 3.8;
- compatibility goal: preserve the established codec/data-structure contract while modernizing implementation and tooling.

A Python-floor increase is reserved for a new major release.

## 1.x stabilization

Before cutting a maintained 1.x release:

1. run the full CI gate on Python 3.8;
2. run the representative newer-Python compatibility gate;
3. confirm Cython and standalone CMake builds;
4. confirm binary serialization round-trip tests;
5. update CHANGELOG.md and README.md;
6. verify pyproject.toml version and requires-python together.

When 1.x reaches a stable release point, create a maintenance branch such as release/1.x if continued Python 3.8 fixes are needed while main moves toward the next major release.

## Future major release

Do not start the next major line merely to obtain newer dependency versions. First finish the 1.x correctness and regression baseline.

When a new Python floor is intentionally selected:

1. create the major-version development branch;
2. update requires-python and CI together;
3. document the final Python-3.8-compatible GVC release;
4. modernize dependency ranges around the new floor;
5. keep format/API breakage explicit and separately tested.
