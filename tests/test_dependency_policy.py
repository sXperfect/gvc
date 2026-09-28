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
