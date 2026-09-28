from __future__ import print_function

import argparse
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CMAKE_SOURCE = PROJECT_ROOT / "library" / "libgvc" / "src"
CMAKE_BUILD = PROJECT_ROOT / "tmp" / "libgvc-build"


def _run(command):
    print("$ " + " ".join(str(part) for part in command), flush=True)
    return subprocess.call(list(command), cwd=str(PROJECT_ROOT))


def run_metadata():
    text = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    version_text = (PROJECT_ROOT / "gvc" / "_version.py").read_text(encoding="utf-8")
    requires = re.search(r'requires-python\s*=\s*"([^"]+)"', text)
    version = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)
    if not requires or requires.group(1) != ">=3.8,<3.13":
        print("unexpected Python support window", file=sys.stderr)
        return 2
    if not version or not version.group(1).startswith("1.0."):
        print("v1.0 development requires a 1.0.x version", file=sys.stderr)
        return 2
    if "tspsolve" in text.lower():
        print("obsolete tspsolve dependency is still declared", file=sys.stderr)
        return 2
    return 0


def run_syntax():
    return _run([
        sys.executable,
        "-m",
        "compileall",
        "-q",
        "gvc",
        "tests",
        "scripts",
    ])


def run_native():
    status = _run([
        sys.executable,
        "-c",
        "import gvc.cquery, gvc.cdebinarize, gvc.data_structures.crc_id",
    ])
    if status:
        return status
    status = _run([
        "cmake",
        "-S",
        str(CMAKE_SOURCE),
        "-B",
        str(CMAKE_BUILD),
        "-DCMAKE_BUILD_TYPE=Release",
    ])
    if status:
        return status
    return _run(["cmake", "--build", str(CMAKE_BUILD), "--parallel", "2"])


def run_test():
    return _run([sys.executable, "-m", "pytest"])


def run_cli():
    return _run([sys.executable, "-m", "gvc", "--help"])


def run_dependencies():
    return _run([
        sys.executable,
        "-c",
        (
            "import sys, numpy, scipy; "
            "print('python', sys.version.split()[0]); "
            "print('numpy', numpy.__version__); "
            "print('scipy', scipy.__version__); "
            "import cyvcf2; print('cyvcf2', cyvcf2.__version__); "
            "import numba; print('numba', numba.__version__)"
        ),
    ])


GATES = {
    "metadata": run_metadata,
    "syntax": run_syntax,
    "native": run_native,
    "test": run_test,
    "cli": run_cli,
    "dependencies": run_dependencies,
}


def run_all():
    for name in ("metadata", "syntax", "native", "test", "cli", "dependencies"):
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
        return run_all()
    return GATES[args.gate]()


if __name__ == "__main__":
    raise SystemExit(main())
