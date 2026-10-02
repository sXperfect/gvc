"""Shared version policy for the GVC 1.0.x maintenance line."""

import re

DEV_RE = re.compile(r"^1\.0\.(\d+)\.dev(\d+)$")
RC_RE = re.compile(r"^1\.0\.(\d+)rc([1-9]\d*)$")
FINAL_RE = re.compile(r"^1\.0\.(\d+)$")


def classify_version(version):
    """Return dev, rc, final, or invalid for a GVC 1.0.x version string."""
    if DEV_RE.fullmatch(version):
        return "dev"
    if RC_RE.fullmatch(version):
        return "rc"
    if FINAL_RE.fullmatch(version):
        return "final"
    return "invalid"


def is_release_line(version):
    return classify_version(version) != "invalid"


def is_taggable(version):
    return classify_version(version) in ("rc", "final")


def expected_tag(version):
    if not is_taggable(version):
        raise ValueError("version is not taggable: {}".format(version))
    return "v" + version


def _parts(version):
    state = classify_version(version)
    if state == "invalid":
        raise ValueError("invalid GVC 1.0.x version: {}".format(version))
    if state == "dev":
        match = DEV_RE.fullmatch(version)
        return state, int(match.group(1)), int(match.group(2))
    if state == "rc":
        match = RC_RE.fullmatch(version)
        return state, int(match.group(1)), int(match.group(2))
    match = FINAL_RE.fullmatch(version)
    return state, int(match.group(1)), None


def validate_transition(current, target):
    """Validate an allowed promotion within one 1.0.x patch release."""
    current_state, current_patch, current_serial = _parts(current)
    target_state, target_patch, target_serial = _parts(target)

    if current_patch != target_patch:
        raise ValueError(
            "version transition must preserve patch number: {} -> {}".format(
                current, target
            )
        )

    if current_state == "dev":
        if target_state != "rc":
            raise ValueError("development versions may only promote to an RC")
        if target_serial is None or target_serial < 1:
            raise ValueError("RC number must be positive")
        return True

    if current_state == "rc":
        if target_state == "rc":
            if target_serial <= current_serial:
                raise ValueError("RC number must increase")
            return True
        if target_state == "final":
            return True
        raise ValueError("RC versions may only promote to a newer RC or final")

    raise ValueError("final versions cannot be promoted in-place")
