#!/usr/bin/env python3
"""Verify that a proposed release tag matches the package version."""

import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def package_version():
    namespace = {}
    exec(
        (ROOT / "gvc" / "_version.py").read_text(encoding="utf-8"),
        namespace,
    )
    return namespace["__version__"]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Check GVC release tag/version consistency."
    )
    parser.add_argument("tag", help="proposed tag, for example v1.0.1rc1")
    args = parser.parse_args(argv)

    version = package_version()
    expected = "v" + version
    if re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", version) is None:
        parser.error(
            "package version must be a 1.0.x RC or final release, got {}".format(
                version
            )
        )
    if args.tag != expected:
        parser.error(
            "tag {} does not match package version {}; expected {}".format(
                args.tag, version, expected
            )
        )

    print("{} matches package version {}".format(args.tag, version))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
