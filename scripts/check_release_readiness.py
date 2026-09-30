#!/usr/bin/env python3
"""Static release-readiness preflight for the GVC 1.0.x line."""

import argparse
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

REQUIRED_PATHS = (
    "docs/audits/release-1.0-readiness.md",
    "docs/audits/python38-v1-compatibility.md",
    "docs/design/multiprocessing.md",
    "benchmarks/run_release.py",
    "benchmarks/compare_release.py",
    "scripts/verify_historical.py",
    "scripts/run_release_validation.py",
    ".github/workflows/release-validation.yml",
    "tests/test_historical_luh_fixture.py",
    "tests/test_jbigkit.py",
    "tests/test_jbigkit_multiprocessing.py",
)

REQUIRED_WORKFLOW_SNIPPETS = (
    "Install JBIG-KIT release dependency",
    "Fetch pinned LUH historical fixture",
    "Verify reviewed Python 3.8 latest stack",
    "Verify Python 3.14 current stack",
)


def check_readiness(allow_dev=True):
    problems = []

    for relative in REQUIRED_PATHS:
        if not (ROOT / relative).is_file():
            problems.append("missing required release file: {}".format(relative))

    version_text = (ROOT / "gvc" / "_version.py").read_text(encoding="utf-8")
    match = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)
    if not match:
        problems.append("could not determine package version")
    else:
        version = match.group(1)
        if not version.startswith("1.0."):
            problems.append("release/1.0 requires a 1.0.x version")
        if not allow_dev:
            rc_or_final = re.fullmatch(r"1\\.0\\.\\d+(?:rc\\d+)?", version)
            if rc_or_final is None:
                problems.append(
                    "release candidate must use 1.0.<patch>rcN or 1.0.<patch>"
                )

    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    if 'requires-python = ">=3.8"' not in pyproject:
        problems.append("release/1.0 must retain Python >=3.8")

    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    for snippet in REQUIRED_WORKFLOW_SNIPPETS:
        if snippet not in workflow:
            problems.append("CI workflow missing gate: {}".format(snippet))

    if "f9af2127a2ff0b87727924860e33f1905fd507cd" not in workflow:
        problems.append("historical LUH fixture is not pinned to reviewed commit")
    if "af45a419e46563906ac51fad0be869291316cea9" not in workflow:
        problems.append("historical LUH fixture blob identity is not checked")

    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Check static prerequisites for the GVC 1.0 release line."
    )
    parser.add_argument(
        "--rc",
        action="store_true",
        help="require an RC/final version instead of allowing .dev",
    )
    args = parser.parse_args(argv)

    problems = check_readiness(allow_dev=not args.rc)
    if problems:
        for problem in problems:
            print("ERROR: " + problem, file=sys.stderr)
        return 1

    print("GVC 1.0 static release-readiness checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
