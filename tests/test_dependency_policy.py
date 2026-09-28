import sys

import numpy as np
import scipy


def _version_tuple(value):
    parts = []
    for token in value.split("."):
        digits = "".join(ch for ch in token if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def test_runtime_dependency_window_matches_python():
    numpy_version = _version_tuple(np.__version__)
    scipy_version = _version_tuple(scipy.__version__)

    if sys.version_info < (3, 9):
        assert numpy_version < (1, 25)
        assert scipy_version < (1, 11)
    else:
        assert numpy_version < (2, 0)
        assert scipy_version < (2, 0)


def test_supported_interpreter_window():
    assert (3, 8) <= sys.version_info[:2] < (3, 13)


def test_optional_dependency_versions_match_python():
    import cyvcf2
    import numba

    cyvcf2_version = _version_tuple(cyvcf2.__version__)
    numba_version = _version_tuple(numba.__version__)

    if sys.version_info < (3, 9):
        assert (0, 31, 4) <= cyvcf2_version < (0, 32)
        assert (0, 57, 1) <= numba_version < (0, 58)
    elif sys.version_info < (3, 10):
        assert (0, 34) <= cyvcf2_version < (0, 35)
        assert (0, 60) <= numba_version < (0, 61)
    else:
        assert (0, 34) <= cyvcf2_version < (0, 35)
        assert (0, 67) <= numba_version < (0, 68)
