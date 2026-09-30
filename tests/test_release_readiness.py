from scripts.check_release_readiness import check_readiness


def test_release_readiness_static_preflight_passes_for_development_tree():
    assert check_readiness(allow_dev=True) == []


def test_release_readiness_rc_mode_rejects_dev_version():
    problems = check_readiness(allow_dev=False)
    assert any("rcN" in problem for problem in problems)
