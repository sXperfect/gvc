import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "check_release_tag",
    ROOT / "scripts" / "check_release_tag.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_package_version_reader_matches_gvc():
    import gvc

    assert MODULE.package_version() == gvc.__version__


@pytest.mark.parametrize(
    "version,tag",
    [
        ("1.0.1rc1", "v1.0.1rc1"),
        ("1.0.1", "v1.0.1"),
    ],
)
def test_release_tag_pattern(version, tag):
    import re

    assert re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", version)
    assert tag == "v" + version


@pytest.mark.parametrize(
    "version",
    ["1.0.1.dev0", "1.1.0rc1", "2.0.0", "1.0.1.post1"],
)
def test_non_release_line_versions_are_rejected(version):
    import re

    assert re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", version) is None
