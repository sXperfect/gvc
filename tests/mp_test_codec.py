"""Small deterministic codec used by multiprocessing integration tests."""

from io import BytesIO

import numpy as np

from gvc.codec.jbig import BIE_HEADER_LEN


def encode(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError("test codec only accepts matrices")
    header = bytearray(BIE_HEADER_LEN)
    header[4:8] = int(matrix.shape[1]).to_bytes(4, "big")
    header[8:12] = int(matrix.shape[0]).to_bytes(4, "big")
    payload = BytesIO()
    np.save(payload, matrix, allow_pickle=False)
    return bytes(header) + payload.getvalue()


def decode(payload):
    payload = bytes(payload)
    if len(payload) < BIE_HEADER_LEN:
        raise ValueError("test codec payload is truncated")
    return np.load(BytesIO(payload[BIE_HEADER_LEN:]), allow_pickle=False)



def install():
    """Install the deterministic test codec inside the current process."""
    from gvc.codec import MAT_CODECS
    from gvc.data_structures.consts import CodecID

    MAT_CODECS[CodecID.JBIG1]["encoder"] = encode
    MAT_CODECS[CodecID.JBIG1]["decoder"] = decode



def signal_child(marker_path):
    """Importable child target for signal-cleanup multiprocessing tests."""
    import os
    import time
    from pathlib import Path

    marker = Path(marker_path)
    marker.write_text(str(os.getpid()))
    while True:
        time.sleep(1)
