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


def test_shared_version_policy_classifies_release_states():
    policy = MODULE._load_version_policy()

    assert policy.classify_version("1.0.1.dev0") == "dev"
    assert policy.classify_version("1.0.1rc1") == "rc"
    assert policy.classify_version("1.0.1") == "final"


@pytest.mark.parametrize(
    "version",
    ["1.1.0rc1", "2.0.0", "1.0.1.post1", "1.0.1a1", "1.0"],
)
def test_shared_version_policy_rejects_invalid_release_line_versions(version):
    policy = MODULE._load_version_policy()

    assert policy.classify_version(version) == "invalid"
    assert not policy.is_release_line(version)
    assert not policy.is_taggable(version)


@pytest.mark.parametrize(
    "version,tag",
    [
        ("1.0.1rc1", "v1.0.1rc1"),
        ("1.0.1rc2", "v1.0.1rc2"),
        ("1.0.1", "v1.0.1"),
        ("1.0.12", "v1.0.12"),
    ],
)
def test_taggable_versions_have_exact_expected_tag(version, tag):
    policy = MODULE._load_version_policy()

    assert policy.is_taggable(version)
    assert policy.expected_tag(version) == tag


@pytest.mark.parametrize("version", ["1.0.1.dev0", "1.0.1.dev3", "1.1.0rc1"])
def test_non_taggable_versions_do_not_have_release_tags(version):
    policy = MODULE._load_version_policy()

    assert not policy.is_taggable(version)
    with pytest.raises(ValueError, match="not taggable"):
        policy.expected_tag(version)
