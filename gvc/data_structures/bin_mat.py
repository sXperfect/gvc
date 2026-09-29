import numpy as np

from .. import utils
from .consts import NCOL_LEN, NROW_LEN


class BinMat:
    """Small binary-matrix framing helper used by legacy GVC structures."""

    def __init__(self, bin_mat):
        matrix = np.asarray(bin_mat)
        if matrix.ndim not in (1, 2):
            raise ValueError("binary matrix must be one- or two-dimensional")
        if matrix.size and np.any((matrix != 0) & (matrix != 1)):
            raise ValueError("binary matrix entries must be 0 or 1")
        self.bin_mat = matrix.astype(np.bool_, copy=False)

    @staticmethod
    def split_data(data):
        data = bytes(data)
        header_len = NROW_LEN + NCOL_LEN
        if len(data) < header_len:
            raise ValueError("binary matrix payload is truncated before dimensions")

        nrows = utils.bytes_to_int(data[:NROW_LEN])
        ncols = utils.bytes_to_int(data[NROW_LEN:header_len])
        payload = data[header_len:]

        expected = nrows if ncols == 0 else nrows * ncols
        if len(payload) != expected:
            raise ValueError(
                "binary matrix payload length mismatch: expected {}, got {}".format(
                    expected, len(payload)
                )
            )
        return nrows, ncols, payload

    def to_bytes(self):
        if self.bin_mat.ndim == 2:
            nrows, ncols = self.bin_mat.shape
        else:
            nrows = self.bin_mat.shape[0]
            ncols = 0

        payload = bytearray()
        payload += utils.int_to_bytes(nrows, NROW_LEN)
        payload += utils.int_to_bytes(ncols, NCOL_LEN)
        payload += self.bin_mat.astype(np.uint8, copy=False).tobytes(order="C")
        return bytes(payload)

    @classmethod
    def from_bytes(cls, data):
        nrows, ncols, payload = cls.split_data(data)
        values = np.frombuffer(payload, dtype=np.uint8)
        if values.size and np.any(values > 1):
            raise ValueError("binary matrix payload contains a non-binary value")

        if ncols == 0:
            matrix = values.reshape(nrows)
        else:
            matrix = values.reshape(nrows, ncols)
        return cls(matrix)
