"""Optional ctypes bridge to the standalone libgvc C helper."""

import ctypes as ct
import os
from pathlib import Path

import numpy as np

_LIBDS = None


def _default_library_path():
    return (
        Path(__file__).resolve().parent.parent
        / "library"
        / "libgvc"
        / "build"
        / "libds.so"
    )


def _load_library():
    global _LIBDS
    if _LIBDS is not None:
        return _LIBDS

    configured = os.environ.get("GVC_LIBDS_PATH")
    path = Path(configured) if configured else _default_library_path()
    if not path.is_file():
        raise FileNotFoundError(
            "libgvc helper not found at {}. Build it with CMake or set "
            "GVC_LIBDS_PATH.".format(path)
        )

    lib = ct.cdll.LoadLibrary(str(path))
    lib.decode_ids_checked.argtypes = [
        ct.POINTER(ct.c_uint8),
        ct.c_size_t,
        ct.POINTER(ct.c_uint16),
        ct.c_uint16,
    ]
    lib.decode_ids_checked.restype = ct.c_int
    _LIBDS = lib
    return lib


def _required_payload_bytes(num_entries):
    if num_entries <= 1:
        return 0
    bits_per_id = (num_entries - 1).bit_length()
    return (num_entries * bits_per_id + 7) // 8


def decode_rowcolids(payload, num_entries):
    if not isinstance(num_entries, int) or isinstance(num_entries, bool):
        raise TypeError("num_entries must be an integer")
    if not 0 <= num_entries <= np.iinfo(np.uint16).max:
        raise ValueError("num_entries is outside the libgvc uint16 range")

    payload = bytes(payload)
    required = _required_payload_bytes(num_entries)
    if len(payload) < required:
        raise ValueError(
            "permutation payload is truncated: expected at least {}, got {}".format(
                required, len(payload)
            )
        )

    recon_ids = np.zeros(num_entries, dtype=np.uint16)
    recon_ids_ptr = recon_ids.ctypes.data_as(ct.POINTER(ct.c_uint16))

    if payload:
        payload_buffer = (ct.c_uint8 * len(payload)).from_buffer_copy(payload)
        payload_ptr = ct.cast(payload_buffer, ct.POINTER(ct.c_uint8))
    else:
        payload_buffer = None
        payload_ptr = ct.POINTER(ct.c_uint8)()

    status = _load_library().decode_ids_checked(
        payload_ptr,
        len(payload),
        recon_ids_ptr,
        num_entries,
    )
    if status == -1:
        raise RuntimeError("libgvc received an invalid output buffer")
    if status == -2:
        raise ValueError("permutation payload is truncated")
    if status == -3:
        raise ValueError("permutation id is outside the valid range")
    if status != 0:
        raise RuntimeError("libgvc decode failed with status {}".format(status))
    return recon_ids
