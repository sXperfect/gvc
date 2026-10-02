from scripts.check_release_readiness import check_readiness
from scripts import run_release_validation


def test_release_readiness_static_preflight_passes_for_development_tree():
    assert check_readiness(allow_dev=True) == []


def test_release_readiness_rc_mode_rejects_dev_version():
    problems = check_readiness(allow_dev=False)
    assert any("rcN" in problem for problem in problems)



def test_release_readiness_requires_manual_only_release_workflow():
    workflow = (
        run_release_validation.ROOT
        / ".github"
        / "workflows"
        / "release-validation.yml"
    ).read_text(encoding="utf-8")

    assert "workflow_dispatch:" in workflow
    assert "\n  push:" not in workflow
    assert "\n  pull_request:" not in workflow
