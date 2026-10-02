# Versioning and Python support

GVC uses release lines to separate interpreter lifecycle changes from codec and
API compatibility.

| GVC line | Minimum Python | Policy |
| --- | ---: | --- |
| 1.0.x | 3.8 | current compatibility line |
| 1.1.x | 3.9 | possible future floor |
| 1.2.x | 3.10 | possible future floor |

Patch releases within a line never raise its Python minimum. A minor release
may raise the minimum after the prior line remains available for users on the
older interpreter. Major versions remain available for actual breaking
GVC API, serialized-format, or codec changes.

For the detailed dependency and CI policy, see
[../development/python-support.md](../development/python-support.md).


## 1.0.x version states

The maintained 1.0 line uses three explicit version states:

- development: `1.0.<patch>.devN`
- release candidate: `1.0.<patch>rcN`
- final: `1.0.<patch>`

Only RC and final versions are taggable. Their Git tags are exact mirrors of the
package version with a leading `v`, for example `v1.0.1rc1` and
`v1.0.1`.

The intended promotion path for a patch release is:

```text
1.0.1.dev0 -> 1.0.1rc1 -> 1.0.1
```

Additional RCs increment the RC number while preserving the same patch number,
for example `1.0.1rc1 -> 1.0.1rc2`. Post releases, alpha/beta releases, and
cross-line tags such as `1.1.0rc1` are not valid states for the 1.0 release
workflow.

The package version in `gvc/_version.py` is the single source of truth. The
wheel metadata, source-distribution root, installed `gvc.__version__`, and
release tag must all agree with that value.
