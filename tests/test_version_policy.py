from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_v1_python_floor_is_explicit():
    metadata = (ROOT / "pyproject.toml").read_text()
    assert 'version = "1.0.0"' in metadata
    assert 'requires-python = ">=3.8"' in metadata


def test_obsolete_tspsolve_is_not_a_dependency():
    metadata = (ROOT / "pyproject.toml").read_text().lower()
    assert "tspsolve" not in metadata
