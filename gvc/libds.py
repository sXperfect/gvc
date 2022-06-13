# Library Data Structure

import os
import io
import ctypes as ct

import numpy as np

lib_path = os.path.join("library", "libgvc", "build", "libds.so")

libds = ct.cdll.LoadLibrary(lib_path)
libds.decode_ids.argtypes = [
    ct.c_void_p,
    ct.c_void_p,
    ct.c_uint16,
]
libds.decode_ids.restype = None

def decode_rowcolids(payload, num_entries):
    recon_ids = np.zeros(num_entries, dtype=np.uint16)
    recon_ids_ptr = recon_ids.ctypes.data_as(
        # ct.POINTER(ct.c_uint16)
        ct.c_void_p
    )

    libds.decode_ids(
        payload,
        recon_ids_ptr,
        num_entries
    )

    return recon_ids