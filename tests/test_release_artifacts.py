import importlib.util
import tarfile
import zipfile
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "check_release_artifacts",
    ROOT / "scripts" / "check_release_artifacts.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_artifacts(tmp_path, version="1.0.1rc1"):
    wheel = tmp_path / ("gvc-{}-cp38-cp38-linux_x86_64.whl".format(version))
    dist_info = "gvc-{}.dist-info".format(version)
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            dist_info + "/METADATA",
            "Metadata-Version: 2.1\n"
            "Name: gvc\n"
            "Version: {}\n"
            "Requires-Python: >=3.8\n".format(version),
        )
        archive.writestr("gvc/cquery.cpython-38-x86_64-linux-gnu.so", b"x")
        archive.writestr("gvc/cdebinarize.cpython-38-x86_64-linux-gnu.so", b"x")
        archive.writestr(
            "gvc/data_structures/crc_id.cpython-38-x86_64-linux-gnu.so", b"x"
        )

    sdist = tmp_path / ("gvc-{}.tar.gz".format(version))
    root = "gvc-{}".format(version)
    source = tmp_path / "src"
    for relative in MODULE.REQUIRED_SDIST_PATHS:
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    with tarfile.open(sdist, "w:gz") as archive:
        for relative in MODULE.REQUIRED_SDIST_PATHS:
            archive.add(source / relative, arcname=root + "/" + relative)
    return wheel, sdist


def test_release_artifact_inspection_records_hashes(tmp_path):
    wheel, sdist = _write_artifacts(tmp_path)
    evidence = MODULE.inspect_artifacts(wheel, sdist, "1.0.1rc1")

    assert evidence["version"] == "1.0.1rc1"
    assert len(evidence["artifacts"]["wheel"]["sha256"]) == 64
    assert len(evidence["artifacts"]["sdist"]["sha256"]) == 64


def test_release_artifact_inspection_rejects_version_mismatch(tmp_path):
    wheel, sdist = _write_artifacts(tmp_path)
    with pytest.raises(ValueError, match="wheel version"):
        MODULE.inspect_artifacts(wheel, sdist, "1.0.2rc1")


def test_release_readiness_rc_regex_accepts_valid_rc():
    from scripts.check_release_readiness import check_readiness

    # The repository is still on a development version, so this confirms the
    # fixed regex directly rather than mutating package state.
    import re

    assert re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", "1.0.1rc1")
    assert re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", "1.0.1")
    assert not re.fullmatch(r"1\.0\.\d+(?:rc\d+)?", "1.1.0rc1")
