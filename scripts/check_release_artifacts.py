#!/usr/bin/env python3
"""Validate GVC release artifacts without importing them."""

import argparse
import email
import hashlib
import json
import re
import tarfile
import zipfile
from pathlib import Path


REQUIRED_WHEEL_MODULES = (
    "gvc/cquery",
    "gvc/cdebinarize",
    "gvc/data_structures/crc_id",
)
REQUIRED_SDIST_PATHS = (
    "pyproject.toml",
    "setup.py",
    "README.md",
    "LICENSE",
    "NOTICE.md",
    "gvc/__init__.py",
    "gvc/_version.py",
    "gvc/cquery.pyx",
    "gvc/cdebinarize.pyx",
    "gvc/data_structures/crc_id.pyx",
)


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _wheel_metadata(archive):
    names = archive.namelist()
    metadata_names = [
        name for name in names if name.endswith(".dist-info/METADATA")
    ]
    wheel_names = [
        name for name in names if name.endswith(".dist-info/WHEEL")
    ]
    entry_point_names = [
        name for name in names if name.endswith(".dist-info/entry_points.txt")
    ]
    if len(metadata_names) != 1:
        raise ValueError("wheel must contain exactly one dist-info/METADATA")
    if len(wheel_names) != 1:
        raise ValueError("wheel must contain exactly one dist-info/WHEEL")
    if len(entry_point_names) != 1:
        raise ValueError("wheel must contain exactly one dist-info/entry_points.txt")
    return (
        email.message_from_bytes(archive.read(metadata_names[0])),
        email.message_from_bytes(archive.read(wheel_names[0])),
        archive.read(entry_point_names[0]).decode("utf-8"),
        names,
    )


def _wheel_native_module_present(names, module):
    prefix = module + "."
    return any(
        name.startswith(prefix)
        and (
            name.endswith(".so")
            or name.endswith(".pyd")
            or name.endswith(".dll")
            or name.endswith(".dylib")
        )
        for name in names
    )


def inspect_artifacts(wheel_path, sdist_path, expected_version):
    evidence = {
        "schema_version": 1,
        "version": expected_version,
        "artifacts": {},
    }

    with zipfile.ZipFile(wheel_path) as archive:
        metadata, wheel_metadata, entry_points, names = _wheel_metadata(archive)
        if metadata.get("Name") != "gvc":
            raise ValueError("wheel metadata Name must be gvc")
        if metadata.get("Version") != expected_version:
            raise ValueError(
                "wheel version {} does not match {}".format(
                    metadata.get("Version"), expected_version
                )
            )
        if metadata.get("Requires-Python") != ">=3.8":
            raise ValueError("wheel must declare Requires-Python >=3.8")
        if wheel_metadata.get("Root-Is-Purelib") != "false":
            raise ValueError("wheel with native extensions must not be purelib")
        if "gvc = gvc.__main__:main" not in entry_points:
            raise ValueError("wheel is missing the gvc console entry point")

        missing_native = [
            module
            for module in REQUIRED_WHEEL_MODULES
            if not _wheel_native_module_present(names, module)
        ]
        if missing_native:
            raise ValueError(
                "wheel is missing native module(s): {}".format(
                    ", ".join(missing_native)
                )
            )

    with tarfile.open(sdist_path, "r:gz") as archive:
        names = archive.getnames()
        roots = {name.split("/", 1)[0] for name in names if "/" in name}
        if len(roots) != 1:
            raise ValueError("sdist must contain one top-level directory")
        root = next(iter(roots))
        missing = [
            relative
            for relative in REQUIRED_SDIST_PATHS
            if "{}/{}".format(root, relative) not in names
        ]
        if missing:
            raise ValueError(
                "sdist is missing required path(s): {}".format(
                    ", ".join(missing)
                )
            )
        expected_root = "gvc-{}".format(expected_version)
        if root != expected_root:
            raise ValueError(
                "sdist root {} does not match {}".format(root, expected_root)
            )

    for kind, path in (("wheel", wheel_path), ("sdist", sdist_path)):
        evidence["artifacts"][kind] = {
            "filename": path.name,
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }

    return evidence


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Inspect GVC wheel/sdist release artifacts."
    )
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--sdist", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--output")
    args = parser.parse_args(argv)

    if re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", args.version) is None:
        parser.error("--version must be a 1.0.x RC or final version")

    wheel = Path(args.wheel).resolve()
    sdist = Path(args.sdist).resolve()
    if not wheel.is_file() or not sdist.is_file():
        parser.error("wheel and sdist must exist")

    try:
        evidence = inspect_artifacts(wheel, sdist, args.version)
    except (OSError, ValueError, zipfile.BadZipFile, tarfile.TarError) as exc:
        parser.error(str(exc))

    rendered = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
        print(output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
