import math

import numpy as np

from ..bitstream import BitIO


def _bits_per_id(num_entries):
    if not isinstance(num_entries, int) or isinstance(num_entries, bool):
        raise TypeError("num_entries must be an integer")
    if num_entries < 0:
        raise ValueError("num_entries must be non-negative")
    if num_entries <= 1:
        return 0
    return int(math.ceil(math.log2(num_entries)))


def _required_bytes(num_entries):
    bits_per_id = _bits_per_id(num_entries)
    return (num_entries * bits_per_id + 7) // 8


def decode_rowcolids(data, num_entries):
    data = np.asarray(data, dtype=np.uint8)
    if data.ndim != 1:
        raise ValueError("permutation payload must be one-dimensional")

    required = _required_bytes(num_entries)
    if data.size != required:
        if data.size < required:
            raise ValueError("permutation payload is truncated")
        raise ValueError("permutation payload contains trailing bytes")

    ids = np.zeros(num_entries, dtype=np.uint16)
    if num_entries <= 1:
        return ids

    decode_rowcolids_loop(data, ids, num_entries)
    if np.any(ids >= num_entries):
        raise ValueError("permutation id is outside the valid range")
    if np.unique(ids).size != num_entries:
        raise ValueError("permutation payload contains duplicate ids")
    return ids


def decode_rowcolids_loop(data, ids, num_entries):
    bits_per_id = _bits_per_id(num_entries)
    if bits_per_id == 0:
        return

    i_byte = 0
    buffer = 0
    buffer_nbits = 0

    for i_entry in range(num_entries):
        while buffer_nbits < bits_per_id:
            buffer <<= 8
            buffer |= int(data[i_byte])
            i_byte += 1
            buffer_nbits += 8

        buffer_nbits -= bits_per_id
        ids[i_entry] = buffer >> buffer_nbits
        mask = (1 << buffer_nbits) - 1
        buffer &= mask


class RowColIds:
    def __init__(self, ids):
        ids = np.asarray(ids)
        if ids.ndim != 1:
            raise ValueError("permutation ids must be one-dimensional")
        if not np.issubdtype(ids.dtype, np.integer):
            raise TypeError("permutation ids must contain integers")
        if ids.size > np.iinfo(np.uint16).max:
            raise ValueError("permutation has too many entries")

        self.ids = ids.astype(np.uint16, copy=False)
        num_entries = self.ids.size
        if num_entries:
            if np.any(self.ids >= num_entries):
                raise ValueError("permutation id is outside the valid range")
            if np.unique(self.ids).size != num_entries:
                raise ValueError("permutation ids must be unique")

    def to_bitio(self):
        num_entries = self.ids.shape[0]
        bits_per_id = _bits_per_id(num_entries)
        ids_bitio = BitIO()

        for curr_id in self.ids:
            ids_bitio.write(int(curr_id), bits_per_id)

        return ids_bitio

    def to_bytes(self):
        return self.to_bitio().to_bytes(align=True)

    @classmethod
    def from_bytes(cls, data, num_entries):
        ids = decode_rowcolids(np.frombuffer(bytes(data), dtype=np.uint8), num_entries)
        return cls(ids)

    @classmethod
    def from_randomaccesshandler(cls, ra_handler, num_ids):
        return cls.from_bytes(ra_handler.read(), num_ids)
