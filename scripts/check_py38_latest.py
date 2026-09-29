#!/usr/bin/env python3
"""Check that Python 3.8 resolves the reviewed dependency ceiling."""

from __future__ import print_function

import importlib.metadata as metadata
import sys

EXPECTED = {
    "numpy": "1.24.4",
    "scipy": "1.10.1",
    "cyvcf2": "0.31.4",
    "numba": "0.58.1",
    "Pillow": "10.4.0",
    "pytest": "8.3.5",
    "Cython": "3.2.9",
    "build": "1.2.2",
    "build": "1.2.2.post1",
}


def main():
    if sys.version_info[:2] != (3, 8):
        print("py38 ceiling check is only meaningful on Python 3.8", file=sys.stderr)
        return 2

    mismatches = []
    for package, expected in EXPECTED.items():
        actual = metadata.version(package)
        print("{:>10}: {}".format(package, actual))
        if actual != expected:
            mismatches.append((package, expected, actual))

    if mismatches:
        print(
            "\nPython 3.8 resolved a dependency stack different from the "
            "reviewed ceiling:",
            file=sys.stderr,
        )
        for package, expected, actual in mismatches:
            print(
                "  {}: reviewed {}, resolved {}".format(package, expected, actual),
                file=sys.stderr,
            )
        print(
            "Review the newer/changed package with the v1.0 compatibility "
            "suite before updating ci/constraints/py38-latest.txt.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
