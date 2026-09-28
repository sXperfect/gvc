from pathlib import Path

import numpy as np

import gvc
import gvc.common


ROOT = Path(__file__).resolve().parents[1]


def test_v1_0_python_floor_is_explicit():
    metadata = (ROOT / "pyproject.toml").read_text()
    assert 'requires-python = ">=3.8"' in metadata
    assert gvc.__version__.startswith("1.0.")


def test_v1_0_dependency_policy_is_modern_and_bounded():
    metadata = (ROOT / "pyproject.toml").read_text()
    expected = [
        '"numpy>=1.24.4,<3"',
        '"scipy>=1.10.1,<2"',
        '"cyvcf2>=0.31.4,<0.32; python_version < \'3.9\'"',
        '"cyvcf2>=0.34.0,<1; python_version >= \'3.9\'"',
        '"numba>=0.58.1,<1"',
        '"Pillow>=10.4.0,<13"',
        '"pytest>=8.3.5,<10"',
        '"Cython>=3.2.9,<4"',
    ]
    for requirement in expected:
        assert requirement in metadata


def test_obsolete_tspsolve_is_not_a_dependency():
    metadata = (ROOT / "pyproject.toml").read_text().lower()
    assert "tspsolve" not in metadata


def test_optional_dependencies_are_not_core_requirements():
    metadata = (ROOT / "pyproject.toml").read_text()
    core = metadata.split("[project.optional-dependencies]", 1)[0].lower()
    assert "cyvcf2" not in core
    assert "numba" not in core
    assert "pillow" not in core


def test_removed_numpy_bool_alias_is_not_used_for_public_dtypes():
    assert gvc.common.PHASING_DTYPE is np.bool_
    assert gvc.common.BIN_DTYPE is np.bool_
