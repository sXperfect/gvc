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
