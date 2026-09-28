from pathlib import Path

import numpy as np

import gvc.common


ROOT = Path(__file__).resolve().parents[1]


def test_v1_0_python_floor_is_explicit():
    metadata = (ROOT / "pyproject.toml").read_text()
    version_file = (ROOT / "gvc" / "_version.py").read_text()

    assert 'requires-python = ">=3.8"' in metadata
    assert '__version__ = "1.0.0.dev0"' in version_file


def test_v1_0_dependency_floors_match_latest_python38_generation():
    metadata = (ROOT / "pyproject.toml").read_text()

    assert '"numpy>=1.24.4"' in metadata
    assert '"scipy>=1.10.1"' in metadata
    assert '"cyvcf2>=0.33.0"' in metadata
    assert '"numba>=0.58.1"' in metadata
    assert '"Pillow>=10.4.0"' in metadata
    assert '"pytest>=8.3.5"' in metadata


def test_obsolete_tspsolve_is_not_a_dependency():
    metadata = (ROOT / "pyproject.toml").read_text().lower()
    assert "tspsolve" not in metadata


def test_removed_numpy_bool_alias_is_not_used_for_public_dtypes():
    assert gvc.common.PHASING_DTYPE is np.bool_
    assert gvc.common.BIN_DTYPE is np.bool_
