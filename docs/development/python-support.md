# Python support policy

## GVC 1.x

GVC 1.x has a minimum supported Python version of **3.8**.

The 1.x line may modernize dependencies, packaging, tests, CI, documentation, and internal implementation as long as those changes do not intentionally raise the Python floor or break the serialized codec contract.

Dependency specifications should use ranges and environment markers so pip can select versions compatible with Python 3.8 rather than freezing the whole project to one historical dependency set.

## Raising the Python floor

A higher minimum Python version is a **major-release boundary** for this repository.

| Release line | Minimum Python |
| --- | ---: |
| 1.x | 3.8 |
| 2.x | to be selected when 2.0 work starts |

When a major line raises the floor:

1. the final compatible 1.x release remains installable for Python 3.8 users;
2. requires-python in the new major line is updated;
3. CI for the new line drops unsupported interpreters;
4. a maintenance branch may be created for critical 1.x fixes;
5. dependency ranges may then be modernized around the new interpreter baseline.

A Python-floor increase must not be hidden inside a patch or minor release.

## CI interpretation

CI should always include the minimum supported interpreter because that is the environment most likely to expose accidental compatibility regressions.

For 1.x, Python 3.8 is therefore a required gate. A newer representative interpreter is also tested to prevent the compatibility line from becoming artificially tied to Python 3.8-era dependency behavior.
