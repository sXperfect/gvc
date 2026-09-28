#!/usr/bin/env python3
"""Local CI dispatcher shared by developers and GitHub Actions."""

from __future__ import print_function

import argparse
import compileall
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"
VERSION_FILE = ROOT / "gvc" / "_version.py"
LIBGVC_BUILD = ROOT / "library" / "libgvc" / "build"


def run(cmd):
    print("$ " + " ".join(str(x) for x in cmd), flush=True)
    return subprocess.call([str(x) for x in cmd], cwd=str(ROOT))


def metadata():
    pyproject = PYPROJECT.read_text(encoding="utf-8")
    version_text = VERSION_FILE.read_text(encoding="utf-8")

    requires = re.search(r'requires-python\s*=\s*"([^"]+)"', pyproject)
    version = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)

    if not requires or requires.group(1) != ">=3.8":
        print("expected requires-python = \">=3.8\"", file=sys.stderr)
        return 1
    if not version or not version.group(1).startswith("1.0."):
        print("GVC 1.0 development must use a 1.0.x version", file=sys.stderr)
        return 1
    return 0


def syntax():
    ok = compileall.compile_dir(str(ROOT / "gvc"), quiet=1)
    ok = compileall.compile_dir(str(ROOT / "tests"), quiet=1) and ok
    ok = compileall.compile_dir(str(ROOT / "scripts"), quiet=1) and ok
    return 0 if ok else 1


def native():
    status = run([sys.executable, "setup.py", "build_ext", "--inplace"])
    if status:
        return status
    status = run(["cmake", "-S", "library/libgvc/src", "-B", str(LIBGVC_BUILD)])
    if status:
        return status
    return run(["cmake", "--build", str(LIBGVC_BUILD), "--parallel"])


def test():
    return run([
        sys.executable,
        "-m",
        "unittest",
        "discover",
        "--start-directory",
        "tests",
        "--verbose",
    ])


def cli():
    return run([sys.executable, "-m", "gvc", "--help"])


GATES = {
    "metadata": metadata,
    "syntax": syntax,
    "native": native,
    "test": test,
    "cli": cli,
}


def all_gates():
    for name in ("metadata", "syntax", "native", "test", "cli"):
        print("\n== {} ==".format(name), flush=True)
        status = GATES[name]()
        if status:
            return status
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gate", choices=list(GATES) + ["all"])
    args = parser.parse_args()
    if args.gate == "all":
        return all_gates()
    return GATES[args.gate]()


if __name__ == "__main__":
    raise SystemExit(main())
