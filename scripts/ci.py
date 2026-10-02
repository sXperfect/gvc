#!/usr/bin/env python3
"""Local CI dispatcher shared by developers and GitHub Actions."""

from __future__ import print_function

import argparse
import compileall
import importlib.metadata as metadata
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from scripts import check_release_artifacts

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"
VERSION_FILE = ROOT / "gvc" / "_version.py"
CMAKE_SOURCE = ROOT / "library" / "libgvc" / "src"
CMAKE_BUILD = ROOT / "tmp" / "libgvc-build"
DIST_DIR = ROOT / "tmp" / "dist"
PACKAGE_VENV_ROOT = ROOT / "tmp" / "package-venvs"

import tarfile
import zipfile

DEPENDENCIES = (
    "numpy",
    "scipy",
    "Cython",
    "pytest",
    "cyvcf2",
    "numba",
    "Pillow",
    "build",
)


def run(cmd):
    print("$ " + " ".join(str(x) for x in cmd), flush=True)
    return subprocess.call([str(x) for x in cmd], cwd=str(ROOT))


def _package_version_for_ci():
    namespace = {}
    exec(VERSION_FILE.read_text(encoding="utf-8"), namespace)
    return namespace["__version__"]


def metadata_gate():
    text = PYPROJECT.read_text(encoding="utf-8")
    version_text = VERSION_FILE.read_text(encoding="utf-8")

    requires = re.search(r'requires-python\s*=\s*"([^"]+)"', text)
    version = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)
    if not requires or requires.group(1) != ">=3.8":
        print('expected requires-python = ">=3.8"', file=sys.stderr)
        return 2
    if not version or not version.group(1).startswith("1.0."):
        print("GVC 1.0 development requires a 1.0.x version", file=sys.stderr)
        return 2

    required = (
        '"numpy>=1.24.4,<3"',
        '"scipy>=1.10.1,<2"',
        '"cyvcf2>=0.31.4,<0.32; python_version < \'3.9\'"',
        '"cyvcf2>=0.34.0,<1; python_version >= \'3.9\'"',
        '"numba>=0.58.1,<1"',
        '"Pillow>=10.4.0,<13"',
        '"pytest>=8.3.5,<10"',
        '"Cython>=3.2.9,<4"',
        '"build>=1.2.2.post1,<2"',
        '"setuptools>=75.3.2,<76; python_version < \'3.9\'"',
        '"wheel>=0.45.1,<1; python_version < \'3.9\'"',
    )
    missing = [item for item in required if item not in text]
    if missing:
        print("missing compatibility policy: " + ", ".join(missing), file=sys.stderr)
        return 2
    if "tspsolve" in text.lower():
        print("obsolete tspsolve dependency is still declared", file=sys.stderr)
        return 2
    return 0


def syntax_gate():
    ok = compileall.compile_dir(str(ROOT / "gvc"), quiet=1)
    ok = compileall.compile_dir(str(ROOT / "tests"), quiet=1) and ok
    ok = compileall.compile_dir(str(ROOT / "scripts"), quiet=1) and ok
    ok = compileall.compile_dir(str(ROOT / "benchmarks"), quiet=1) and ok
    return 0 if ok else 1


def native_gate():
    # The environment is prepared with `pip install -e .`, which builds the
    # Cython extensions through the declared PEP 517 backend. Verify those
    # installed/in-place extension modules directly instead of invoking the
    # deprecated setup.py command path.
    status = run([
        sys.executable,
        "-c",
        "import gvc.cquery, gvc.cdebinarize, gvc.data_structures.crc_id",
    ])
    if status:
        return status

    if CMAKE_BUILD.exists():
        shutil.rmtree(str(CMAKE_BUILD))
    status = run([
        "cmake",
        "-S",
        str(CMAKE_SOURCE),
        "-B",
        str(CMAKE_BUILD),
        "-DCMAKE_BUILD_TYPE=Release",
    ])
    if status:
        return status
    return run(["cmake", "--build", str(CMAKE_BUILD), "--parallel", "2"])


def test_gate():
    return run([sys.executable, "-m", "pytest"])


def optional_gate():
    code = (
        "import cyvcf2, numba, PIL; "
        "print('cyvcf2', cyvcf2.__version__); "
        "print('numba', numba.__version__); "
        "print('Pillow', PIL.__version__)"
    )
    return run([sys.executable, "-c", code])



def _venv_python(venv_dir):
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _run_external(cmd, cwd, env=None):
    print("$ " + " ".join(str(x) for x in cmd), flush=True)
    return subprocess.call(
        [str(x) for x in cmd],
        cwd=str(cwd),
        env=env,
    )


def _install_and_smoke_artifact(artifact, name):
    venv_dir = PACKAGE_VENV_ROOT / name
    if venv_dir.exists():
        shutil.rmtree(str(venv_dir))
    status = run([sys.executable, "-m", "venv", str(venv_dir)])
    if status:
        return status

    py = _venv_python(venv_dir)
    env = dict(os.environ)
    if sys.version_info[:2] == (3, 8):
        env["PIP_CONSTRAINT"] = str(ROOT / "ci" / "constraints" / "py38-latest.txt")

    pip_spec = "pip<25.1" if sys.version_info[:2] == (3, 8) else "pip"
    status = _run_external(
        [py, "-m", "pip", "install", "--upgrade", pip_spec],
        cwd=venv_dir,
        env=env,
    )
    if status:
        return status

    status = _run_external(
        [py, "-m", "pip", "install", str(artifact)],
        cwd=venv_dir,
        env=env,
    )
    if status:
        return status

    outside = Path(tempfile.mkdtemp(prefix="gvc-package-smoke-"))
    code = (
        "from pathlib import Path; "
        "import importlib.metadata as md; "
        "import gvc, gvc.cquery, gvc.cdebinarize, gvc.data_structures.crc_id; "
        "origin = Path(gvc.__file__).resolve(); "
        "print('gvc origin:', origin); "
        "print('gvc version:', md.version('gvc')); "
        "assert md.version('gvc') == gvc.__version__; "
        "source = Path(" + repr(str(ROOT / "gvc")) + ").resolve(); "
        "assert source != origin and source not in origin.parents"
    )
    status = _run_external([py, "-c", code], cwd=outside, env=env)
    if status:
        return status
    return _run_external([py, "-m", "gvc", "--help"], cwd=outside, env=env)


def packaging_gate():
    if DIST_DIR.exists():
        shutil.rmtree(str(DIST_DIR))
    if PACKAGE_VENV_ROOT.exists():
        shutil.rmtree(str(PACKAGE_VENV_ROOT))
    DIST_DIR.mkdir(parents=True, exist_ok=True)
    PACKAGE_VENV_ROOT.mkdir(parents=True, exist_ok=True)

    status = run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--sdist",
            "--outdir",
            str(DIST_DIR),
        ]
    )
    if status:
        return status

    wheels = sorted(DIST_DIR.glob("gvc-*.whl"))
    sdists = sorted(DIST_DIR.glob("gvc-*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        print(
            "expected exactly one wheel and one sdist, got {} wheel(s) and {} sdist(s)".format(
                len(wheels), len(sdists)
            ),
            file=sys.stderr,
        )
        return 2

    try:
        check_release_artifacts.inspect_artifacts(
            wheels[0],
            sdists[0],
            _package_version_for_ci(),
        )
    except (OSError, ValueError, zipfile.BadZipFile, tarfile.TarError) as exc:
        print("artifact inspection failed: {}".format(exc), file=sys.stderr)
        return 2

    status = _install_and_smoke_artifact(wheels[0], "wheel")
    if status:
        return status
    return _install_and_smoke_artifact(sdists[0], "sdist")

def cli_gate():
    return run([sys.executable, "-m", "gvc", "--help"])


def dependency_gate():
    print("Python {}".format(sys.version.replace("\n", " ")))
    for package in DEPENDENCIES:
        try:
            value = metadata.version(package)
        except metadata.PackageNotFoundError:
            value = "<not installed>"
        print("{:>10}: {}".format(package, value))
    return 0


def release_readiness_gate():
    return run([sys.executable, "scripts/check_release_readiness.py"])


GATES = {
    "metadata": metadata_gate,
    "syntax": syntax_gate,
    "native": native_gate,
    "test": test_gate,
    "optional": optional_gate,
    "packaging": packaging_gate,
    "cli": cli_gate,
    "deps": dependency_gate,
    "dependencies": dependency_gate,
    "release": release_readiness_gate,
}


def all_gates():
    for name in ("metadata", "syntax", "native", "test", "cli", "packaging", "release", "optional", "deps"):
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
