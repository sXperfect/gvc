import sys

import numpy as np
import pytest
import scipy


def _version_tuple(value):
    parts = []
    for token in value.split("."):
        digits = "".join(ch for ch in token if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def test_supported_interpreter_floor():
    assert sys.version_info[:2] >= (3, 8)


def test_core_dependency_versions_stay_inside_v1_compatibility():
    numpy_version = _version_tuple(np.__version__)
    scipy_version = _version_tuple(scipy.__version__)

    assert (1, 24, 4) <= numpy_version < (3, 0)
    assert (1, 10, 1) <= scipy_version < (2, 0)


def test_optional_dependency_versions_when_installed():
    cyvcf2 = pytest.importorskip("cyvcf2")
    numba = pytest.importorskip("numba")
    pillow = pytest.importorskip("PIL")

    assert (0, 33, 0) <= _version_tuple(cyvcf2.__version__) < (1, 0)
    assert (0, 58, 1) <= _version_tuple(numba.__version__) < (1, 0)
    assert (10, 4, 0) <= _version_tuple(pillow.__version__) < (13, 0)
