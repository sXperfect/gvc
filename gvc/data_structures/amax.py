import io

import numpy as np

from . import consts
from ..bitstream import BitIO, BitstreamReader


class VectorAMax:
    def __init__(self, vector):
        vector = np.asarray(vector)
        if vector.ndim != 1:
            raise ValueError("AMax vector must be one-dimensional")
        if not np.issubdtype(vector.dtype, np.integer):
            raise TypeError("AMax vector must contain integers")
        if vector.size and np.any(vector < 1):
            raise ValueError("AMax vector entries must be positive and non-zero")
        if vector.size >= (1 << (8 * consts.AMAX_NUM_ENTRIES_LEN)):
            raise ValueError("AMax vector has too many entries")
        self.vector = vector

    @classmethod
    def from_bytes(cls, data):
        data = bytes(data)
        istream = BitstreamReader(io.BytesIO(data))

        num_entries = istream.read_bytes(
            consts.AMAX_NUM_ENTRIES_LEN, ret_int=True
        )
        bits_per_entry = istream.read_bytes(
            consts.AMAX_BITS_PER_ENTRY_LEN, ret_int=True
        )

        header_len = (
            consts.AMAX_NUM_ENTRIES_LEN
            + consts.AMAX_BITS_PER_ENTRY_LEN
        )
        remaining_bits = (len(data) - header_len) * 8
        # Every entry requires at least its presence flag. Reject impossible
        # counts before allocating the output vector so malformed input cannot
        # trigger an attacker-controlled large allocation.
        if num_entries > remaining_bits:
            raise ValueError("AMax payload entry count exceeds available data")

        vector = np.ones(num_entries, dtype=np.uint64)
        for i in range(num_entries):
            flag = istream.read_bits(consts.AMAX_FLAG_BITLEN)
            if flag:
                value = istream.read_bits(bits_per_entry)
                vector[i] = value + 2

        istream.align_to_byte()
        if istream.tell() != len(data):
            raise ValueError("AMax payload contains trailing bytes")

        return cls(vector)

    def to_bitio(self):
        data_bitio = BitIO()
        num_entries = self.vector.shape[0]

        if num_entries:
            max_val = int(np.max(self.vector))
            if max_val != 1:
                offset_max = max_val - 2
                bits_per_entry = int(np.ceil(np.log2(offset_max + 1)))
            else:
                bits_per_entry = 0
        else:
            bits_per_entry = 0

        if bits_per_entry >= (1 << (8 * consts.AMAX_BITS_PER_ENTRY_LEN)):
            raise ValueError("AMax entry width is outside the serialized range")

        data_bitio.write(num_entries, consts.AMAX_NUM_ENTRIES_LEN * 8)
        data_bitio.write(bits_per_entry, consts.AMAX_BITS_PER_ENTRY_LEN * 8)

        for value in self.vector:
            value = int(value)
            flag = value != 1
            data_bitio.write(int(flag), consts.AMAX_FLAG_BITLEN)
            if flag:
                data_bitio.write(value - 2, bits_per_entry)

        return data_bitio
