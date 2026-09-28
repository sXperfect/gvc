#!/usr/bin/env python3
"""Local CI dispatcher shared by developers and GitHub Actions."""

from __future__ import print_function

import argparse
import compileall
import importlib.metadata as metadata
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"
VERSION_FILE = ROOT / "gvc" / "_version.py"
CMAKE_SOURCE = ROOT / "library" / "libgvc" / "src"
CMAKE_BUILD = ROOT / "tmp" / "libgvc-build"

DEPENDENCIES = (
    "numpy",
    "scipy",
    "Cython",
    "pytest",
    "cyvcf2",
    "numba",
    "Pillow",
)


def run(cmd):
    print("$ " + " ".join(str(x) for x in cmd), flush=True)
    return subprocess.call([str(x) for x in cmd], cwd=str(ROOT))


def metadata_gate():
    pyproject = PYPROJECT.read_text(encoding="utf-8")
    version_text = VERSION_FILE.read_text(encoding="utf-8")

    requires = re.search(r'requires-python\s*=\s*"([^"]+)"', pyproject)
    version = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)

    if not requires or requires.group(1) != ">=3.8":
        print('expected requires-python = ">=3.8"', file=sys.stderr)
        return 1
    if not version or not version.group(1).startswith("1.0."):
        print("GVC 1.0 development must use a 1.0.x version", file=sys.stderr)
        return 1

    required_fragments = (
        '"numpy>=1.24.4"',
        '"scipy>=1.10.1"',
        '"cyvcf2>=0.33.0"',
        '"numba>=0.58.1"',
        '"Pillow>=10.4.0"',
        '"pytest>=8.3.5"',
        '"Cython>=3.2.9"',
    )
    missing = [fragment for fragment in required_fragments if fragment not in pyproject]
    if missing:
        print("missing dependency policy: " + ", ".join(missing), file=sys.stderr)
        return 1
    if "tspsolve" in pyproject.lower():
        print("tspsolve must not return as a dependency", file=sys.stderr)
        return 1
    return 0


def syntax_gate():
    ok = compileall.compile_dir(str(ROOT / "gvc"), quiet=1)
    ok = compileall.compile_dir(str(ROOT / "tests"), quiet=1) and ok
    ok = compileall.compile_dir(str(ROOT / "scripts"), quiet=1) and ok
    return 0 if ok else 1


def native_gate():
    status = run([sys.executable, "setup.py", "build_ext", "--inplace"])
    if status:
        return status

    if CMAKE_BUILD.exists():
        shutil.rmtree(str(CMAKE_BUILD))

    status = run(
        [
            "cmake",
            "-S",
            str(CMAKE_SOURCE),
            "-B",
            str(CMAKE_BUILD),
            "-DCMAKE_BUILD_TYPE=Release",
        ]
    )
    if status:
        return status
    return run(["cmake", "--build", str(CMAKE_BUILD), "--parallel", "2"])


def test_gate():
    return run([sys.executable, "-m", "pytest"])


def optional_gate():
    code = (
        "import cyvcf2\n"
        "import numba\n"
        "import PIL\n"
        "print('optional dependency imports OK')\n"
    )
    return run([sys.executable, "-c", code])


def cli_gate():
    return run([sys.executable, "-m", "gvc", "--help"])


def dependency_gate():
    print("Python {}".format(sys.version.replace("\n", " ")))
    for package in DEPENDENCIES:
        try:
            version = metadata.version(package)
        except metadata.PackageNotFoundError:
            version = "<not installed>"
        print("{:>10}: {}".format(package, version))
    return 0


GATES = {
    "metadata": metadata_gate,
    "syntax": syntax_gate,
    "native": native_gate,
    "test": test_gate,
    "optional": optional_gate,
    "cli": cli_gate,
    "deps": dependency_gate,
}


def all_gates():
    for name in ("metadata", "syntax", "native", "test", "cli", "deps"):
        print("\n== {} ==".format(name), flush=True)
        status = GATES[name]()
        if status:
            return status
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="GVC local CI dispatcher")
    parser.add_argument("gate", choices=list(GATES) + ["all"])
    args = parser.parse_args(argv)
    if args.gate == "all":
        return all_gates()
    return GATES[args.gate]()


if __name__ == "__main__":
    raise SystemExit(main())
