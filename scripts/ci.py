from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CMAKE_SOURCE = PROJECT_ROOT / "library" / "libgvc" / "src"
CMAKE_BUILD = PROJECT_ROOT / "tmp" / "libgvc-build"


def _run(command):
    print("$ " + " ".join(str(part) for part in command), flush=True)
    result = subprocess.run(list(command), cwd=str(PROJECT_ROOT))
    return result.returncode


def run_metadata():
    text = (PROJECT_ROOT / "pyproject.toml").read_text()
    required = [
        'version = "1.0.0"',
        'requires-python = ">=3.8"',
    ]
    missing = [entry for entry in required if entry not in text]
    if missing:
        print("missing required metadata: " + ", ".join(missing), file=sys.stderr)
        return 2
    return 0


def run_syntax():
    return _run(
        [
            sys.executable,
            "-m",
            "compileall",
            "-q",
            "gvc",
            "tests",
            "scripts",
        ]
    )


def run_native():
    status = _run(
        [
            sys.executable,
            "-c",
            "import gvc.cquery, gvc.cdebinarize, gvc.data_structures.crc_id",
        ]
    )
    if status != 0:
        return status

    if CMAKE_BUILD.exists():
        shutil.rmtree(str(CMAKE_BUILD))

    status = _run(
        [
            "cmake",
            "-S",
            str(CMAKE_SOURCE),
            "-B",
            str(CMAKE_BUILD),
            "-DCMAKE_BUILD_TYPE=Release",
        ]
    )
    if status != 0:
        return status

    return _run(["cmake", "--build", str(CMAKE_BUILD), "--parallel", "2"])


def run_test():
    return _run([sys.executable, "-m", "pytest"])


GATES = {
    "metadata": run_metadata,
    "syntax": run_syntax,
    "native": run_native,
    "test": run_test,
}


def run_all():
    for name, runner in GATES.items():
        print("\n== {} ==".format(name), flush=True)
        status = runner()
        if status != 0:
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
