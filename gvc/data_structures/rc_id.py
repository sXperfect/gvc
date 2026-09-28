import math

import numpy as np

from ..bitstream import BitIO


def _bits_per_id(num_entries):
    if num_entries < 0:
        raise ValueError("num_entries must be non-negative")
    if num_entries <= 1:
        return 0
    return int(math.ceil(math.log2(num_entries)))


def decode_rowcolids(data, num_entries):
    ids = np.zeros(num_entries, dtype=np.uint16)
    if num_entries == 0:
        return ids
    decode_rowcolids_loop(data, ids, num_entries)
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
            if i_byte >= len(data):
                raise ValueError("permutation payload is truncated")
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
        self.ids = np.asarray(ids, dtype=np.uint16)
        if self.ids.ndim != 1:
            raise ValueError("permutation ids must be one-dimensional")

    def to_bitio(self):
        num_entries = self.ids.shape[0]
        bits_per_id = _bits_per_id(num_entries)

        ids_bitio = BitIO()
        if bits_per_id == 0:
            return ids_bitio

        if np.any(self.ids >= num_entries):
            raise ValueError("permutation id is outside the valid range")

        for curr_id in self.ids:
            ids_bitio.write(int(curr_id), bits_per_id)

        return ids_bitio

    def to_bytes(self):
        return self.to_bitio().to_bytes(align=True)

    @classmethod
    def from_bytes(cls, data, num_entries):
        ids = decode_rowcolids(np.frombuffer(data, dtype=np.uint8), num_entries)
        return cls(ids)

    @classmethod
    def from_randomaccesshandler(cls, ra_handler, num_ids):
        return cls.from_bytes(ra_handler.read(), num_ids)
