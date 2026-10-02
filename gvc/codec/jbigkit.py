"""Maintained JBIG-KIT subprocess integration for the historical GVC codec."""

import os
import shutil
import subprocess
import tempfile

import numpy as np

from .jbig import BIE_HEADER_LEN, get_shape


DEFAULT_ENCODER = "pbmtojbg85"
DEFAULT_DECODER = "jbgtopbm85"
ENCODER_ENV = "GVC_JBIG_ENCODER"
DECODER_ENV = "GVC_JBIG_DECODER"
TIMEOUT_ENV = "GVC_JBIG_TIMEOUT"
DEFAULT_TIMEOUT_SECONDS = 60.0

JBG_TPBON = 0x08


class JBIGKitError(RuntimeError):
    """Raised when an external JBIG-KIT command fails."""


def _pillow_image():
    try:
        from PIL import Image
    except ImportError as exc:
        raise ImportError(
            "JBIG-KIT integration requires Pillow; install gvc[jbig] or gvc[all]"
        ) from exc
    Image.MAX_IMAGE_PIXELS = None
    return Image


def _resolve_executable(explicit, env_name, default_name):
    candidate = explicit or os.environ.get(env_name)
    if candidate:
        candidate = os.path.abspath(os.path.expanduser(candidate))
        if not os.path.isfile(candidate):
            raise FileNotFoundError(
                "{} executable does not exist: {}".format(env_name, candidate)
            )
        if not os.access(candidate, os.X_OK):
            raise PermissionError(
                "{} executable is not executable: {}".format(env_name, candidate)
            )
        return candidate

    discovered = shutil.which(default_name)
    if discovered is None:
        raise FileNotFoundError(
            "JBIG-KIT executable {!r} was not found on PATH; install jbigkit-bin "
            "or set {}".format(default_name, env_name)
        )
    return discovered


def _timeout(value):
    if value is None:
        raw = os.environ.get(TIMEOUT_ENV)
        value = DEFAULT_TIMEOUT_SECONDS if raw is None else float(raw)
    value = float(value)
    if value <= 0:
        raise ValueError("JBIG timeout must be positive")
    return value


def _run(args, timeout):
    try:
        completed = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise JBIGKitError(
            "JBIG-KIT command timed out after {:.3g}s: {}".format(
                timeout, args[0]
            )
        ) from exc
    except OSError as exc:
        raise JBIGKitError(
            "failed to execute JBIG-KIT command {}: {}".format(args[0], exc)
        ) from exc

    if completed.returncode != 0:
        stderr = completed.stderr.decode("utf-8", errors="replace").strip()
        raise JBIGKitError(
            "JBIG-KIT command failed with exit status {}: {}{}".format(
                completed.returncode,
                args[0],
                ": " + stderr if stderr else "",
            )
        )
    return completed


def _validate_matrix(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError("JBIG matrix must be two-dimensional")
    if matrix.shape[0] <= 0 or matrix.shape[1] <= 0:
        raise ValueError("JBIG matrix dimensions must be positive")
    if matrix.size and np.any((matrix != 0) & (matrix != 1)):
        raise ValueError("JBIG matrix must contain only binary values")
    return matrix.astype(np.bool_, copy=False)


def available(encoder_path=None, decoder_path=None):
    """Return whether both historical T.85 JBIG-KIT executables are usable."""
    try:
        _resolve_executable(encoder_path, ENCODER_ENV, DEFAULT_ENCODER)
        _resolve_executable(decoder_path, DECODER_ENV, DEFAULT_DECODER)
    except (FileNotFoundError, PermissionError):
        return False
    return True


def encode(matrix, encoder_path=None, timeout=None):
    """Encode a binary 2-D matrix using JBIG-KIT pbmtojbg85."""
    matrix = _validate_matrix(matrix)
    encoder = _resolve_executable(
        encoder_path, ENCODER_ENV, DEFAULT_ENCODER
    )
    timeout = _timeout(timeout)
    Image = _pillow_image()

    with tempfile.TemporaryDirectory(prefix="gvc-pbm2jbg-") as directory:
        pbm_path = os.path.join(directory, "matrix.pbm")
        jbg_path = os.path.join(directory, "matrix.jbg")

        image = Image.fromarray(matrix)
        image.save(pbm_path, format="PPM")

        args = [
            encoder,
            "-p",
            str(JBG_TPBON),
            "-s",
            str((1 << 32) - 1),
            pbm_path,
            jbg_path,
        ]
        _run(args, timeout)

        try:
            with open(jbg_path, "rb") as handle:
                payload = handle.read()
        except OSError as exc:
            raise JBIGKitError(
                "JBIG-KIT encoder did not produce a readable output"
            ) from exc

    if len(payload) < BIE_HEADER_LEN:
        raise JBIGKitError(
            "JBIG-KIT encoder produced a truncated payload ({} bytes)".format(
                len(payload)
            )
        )

    encoded_shape = get_shape(payload[:BIE_HEADER_LEN])
    if encoded_shape != matrix.shape:
        raise JBIGKitError(
            "JBIG-KIT encoded shape {} does not match input {}".format(
                encoded_shape, matrix.shape
            )
        )
    return payload


def decode(payload, decoder_path=None, timeout=None):
    """Decode a JBIG1/T.85 payload using JBIG-KIT jbgtopbm85."""
    if not isinstance(payload, (bytes, bytearray, memoryview)):
        raise TypeError("JBIG payload must be bytes-like")
    payload = bytes(payload)
    if len(payload) < BIE_HEADER_LEN:
        raise ValueError(
            "JBIG payload is truncated: expected at least {} bytes".format(
                BIE_HEADER_LEN
            )
        )

    expected_shape = get_shape(payload[:BIE_HEADER_LEN])
    decoder = _resolve_executable(
        decoder_path, DECODER_ENV, DEFAULT_DECODER
    )
    timeout = _timeout(timeout)
    Image = _pillow_image()

    with tempfile.TemporaryDirectory(prefix="gvc-jbg2pbm-") as directory:
        jbg_path = os.path.join(directory, "matrix.jbg")
        pbm_path = os.path.join(directory, "matrix.pbm")

        with open(jbg_path, "wb") as handle:
            handle.write(payload)

        args = [
            decoder,
            "-x",
            str(expected_shape[1]),
            "-B",
            str(1 << 30),
            jbg_path,
            pbm_path,
        ]
        _run(args, timeout)

        try:
            with Image.open(pbm_path) as image:
                matrix = np.asarray(image).copy()
        except (OSError, ValueError) as exc:
            raise JBIGKitError(
                "JBIG-KIT decoder did not produce a readable PBM image"
            ) from exc

    if matrix.shape != expected_shape:
        raise JBIGKitError(
            "JBIG-KIT decoded shape {} does not match payload header {}".format(
                matrix.shape, expected_shape
            )
        )
    return matrix.astype(np.bool_, copy=False)


def install():
    """Register this integration as GVC's JBIG1 codec in the current process."""
    from . import MAT_CODECS
    from ..data_structures.consts import CodecID

    MAT_CODECS[CodecID.JBIG1]["encoder"] = encode
    MAT_CODECS[CodecID.JBIG1]["decoder"] = decode
