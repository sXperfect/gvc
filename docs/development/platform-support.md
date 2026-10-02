# Platform support for GVC 1.0.x

## Validated platform

The GVC 1.0.x release line is validated on **Linux**.

Automated CI, native-extension builds, standalone libgvc checks, historical
compatibility gates, multiprocessing reliability tests, packaging validation,
and controlled release-validation workflows currently execute on Linux.

A 1.0.x release may therefore claim Linux as its tested/validated platform.

## macOS and Windows

macOS and Windows are currently **best-effort** rather than release-gated
platforms.

The codebase contains portability provisions such as Windows compiler flags and
spawn-compatible multiprocessing support, but those provisions do not by
themselves constitute release validation.

Until dedicated CI or equivalent retained release evidence exists, the project
must not claim that macOS or Windows have the same validation level as Linux.

## Promotion to validated status

A platform may be promoted to validated status only after the release process
covers, at minimum:

- wheel/sdist build and isolated installation;
- native Cython extension import;
- standalone native helper build where applicable;
- complete maintained test suite;
- multiprocessing lifecycle behavior appropriate to the platform;
- CLI smoke tests;
- v1 golden-format compatibility;
- a documented external-JBIG strategy or an explicit platform limitation.

Platform support statements in README/release notes must reflect the actual
validated matrix rather than inferred portability from source code.
