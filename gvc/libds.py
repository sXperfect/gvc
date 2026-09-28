"""Optional ctypes bridge to the standalone libgvc C helper."""

import ctypes as ct
import os
from pathlib import Path

import numpy as np

_LIBDS = None


def _default_library_path():
    return Path(__file__).resolve().parent.parent / "library" / "libgvc" / "build" / "libds.so"


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
    lib.decode_ids.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_uint16]
    lib.decode_ids.restype = None
    _LIBDS = lib
    return lib


def decode_rowcolids(payload, num_entries):
    libds = _load_library()
    recon_ids = np.zeros(num_entries, dtype=np.uint16)
    recon_ids_ptr = recon_ids.ctypes.data_as(ct.c_void_p)
    libds.decode_ids(payload, recon_ids_ptr, num_entries)
    return recon_ids
