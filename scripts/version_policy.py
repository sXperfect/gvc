"""Shared version policy for the GVC 1.0.x maintenance line."""

import re

DEV_RE = re.compile(r"^1\.0\.(\d+)\.dev(\d+)$")
RC_RE = re.compile(r"^1\.0\.(\d+)rc(\d+)$")
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
