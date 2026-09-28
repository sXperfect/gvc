import builtins

import numpy as np
import pytest

from gvc import reader
from gvc.dist import hamming_rl_dist


def test_vcf_dependency_is_lazy_and_has_actionable_error(monkeypatch):
    original_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name == "cyvcf2":
            raise ImportError("blocked for test")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked_import)

    with pytest.raises(ImportError, match=r"gvc\[vcf\]"):
        reader._vcf_class()


def test_distance_function_remains_correct_with_or_without_numba():
    left = np.array([0, 0, 1, 1, 0, 1], dtype=np.uint8)
    right = np.array([0, 1, 1, 0, 0, 0], dtype=np.uint8)
    assert hamming_rl_dist(left, right) == 3


@pytest.mark.optional
def test_optional_dependency_stack_imports_when_installed():
    pytest.importorskip("cyvcf2")
    pytest.importorskip("numba")
    pytest.importorskip("PIL")
