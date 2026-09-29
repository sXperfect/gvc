import os
from pathlib import Path

import numpy as np
import pytest

from gvc.codec import jbigkit
from gvc.codec.jbig import BIE_HEADER_LEN, get_shape


def _write_executable(path, body):
    path.write_text(body)
    path.chmod(0o755)
    return str(path)


def test_jbigkit_validates_binary_matrix_and_timeout(tmp_path):
    with pytest.raises(ValueError, match="two-dimensional"):
        jbigkit.encode(np.array([0, 1], dtype=np.uint8), encoder_path="/bin/true")
    with pytest.raises(ValueError, match="binary"):
        jbigkit.encode(
            np.array([[0, 2]], dtype=np.uint8),
            encoder_path="/bin/true",
        )
    with pytest.raises(ValueError, match="positive"):
        jbigkit._timeout(0)


def test_jbigkit_reports_missing_executable(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        jbigkit.encode(
            np.array([[0, 1], [1, 0]], dtype=np.uint8),
            encoder_path=str(tmp_path / "missing"),
        )


def test_jbigkit_reports_subprocess_failure(tmp_path):
    failing = _write_executable(
        tmp_path / "fail.py",
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "sys.stderr.write('synthetic jbig failure\\n')\n"
        "raise SystemExit(7)\n",
    )
    with pytest.raises(jbigkit.JBIGKitError, match="exit status 7"):
        jbigkit.encode(
            np.array([[0, 1], [1, 0]], dtype=np.uint8),
            encoder_path=failing,
            timeout=2,
        )


def test_jbigkit_reports_timeout(tmp_path):
    sleeper = _write_executable(
        tmp_path / "sleep.py",
        "#!/usr/bin/env python3\n"
        "import time\n"
        "time.sleep(5)\n",
    )
    with pytest.raises(jbigkit.JBIGKitError, match="timed out"):
        jbigkit.encode(
            np.array([[0, 1], [1, 0]], dtype=np.uint8),
            encoder_path=sleeper,
            timeout=0.05,
        )


def test_real_jbigkit_roundtrip_when_required():
    if os.environ.get("GVC_REQUIRE_JBIGKIT") != "1":
        pytest.skip("real JBIG-KIT integration is release-gate only")

    if not jbigkit.available():
        pytest.fail(
            "GVC_REQUIRE_JBIGKIT=1 but pbmtojbg85/jbgtopbm85 are unavailable"
        )

    rng = np.random.default_rng(4811)
    matrix = rng.integers(0, 2, size=(17, 29), dtype=np.uint8).astype(bool)

    payload = jbigkit.encode(matrix, timeout=10)
    assert len(payload) >= BIE_HEADER_LEN
    assert get_shape(payload[:BIE_HEADER_LEN]) == matrix.shape

    restored = jbigkit.decode(payload, timeout=10)
    assert restored.dtype == np.bool_
    np.testing.assert_array_equal(restored, matrix)


def test_real_jbigkit_rejects_truncated_payload_when_required():
    if os.environ.get("GVC_REQUIRE_JBIGKIT") != "1":
        pytest.skip("real JBIG-KIT integration is release-gate only")

    with pytest.raises(ValueError, match="truncated"):
        jbigkit.decode(b"\x00" * (BIE_HEADER_LEN - 1))
