from pathlib import Path

import gvc


ROOT = Path(__file__).resolve().parents[1]


def test_v1_python_window_is_explicit():
    metadata = (ROOT / "pyproject.toml").read_text()
    assert 'requires-python = ">=3.8,<3.13"' in metadata
    assert gvc.__version__.startswith("1.0.")


def test_obsolete_tspsolve_is_not_a_dependency():
    metadata = (ROOT / "pyproject.toml").read_text().lower()
    assert "tspsolve" not in metadata


def test_heavy_optional_dependencies_are_not_core():
    metadata = (ROOT / "pyproject.toml").read_text()
    core = metadata.split("[project.optional-dependencies]", 1)[0].lower()
    assert "cyvcf2" not in core
    assert "numba" not in core
    assert "pillow" not in core
