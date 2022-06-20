import io
import logging as log
import typing as t
import numpy as np
from . import consts
from ..bitstream import BitstreamReader, BitIO

class VectorAMax(object):
    def __init__(self, vector):

        assert vector.ndim == 1
        assert np.all(vector != 0)

        self.vector = vector

    @classmethod
    def from_bytes(
        cls, 
        data:bytes
    ):
        istream = BitstreamReader(io.BytesIO(data))

        num_entries = istream.read_bytes(consts.AMAX_NUM_ENTRIES_LEN, ret_int=True)
        bits_per_entry = istream.read_bytes(consts.AMAX_BITS_PER_ENTRY_LEN, ret_int=True)

        vector = np.ones((num_entries), dtype=int)
        for i in range(num_entries):
            flag = istream.read_bits(1)

            if flag:
                value = istream.read_bits(bits_per_entry)
                vector[i] = value+2

        return cls(vector)

    def to_bitio(self):
        data_bitio = BitIO()

        num_entries = self.vector.shape[0]
        max_val = np.max(self.vector)

        if max_val != 1:
            offseted_max_val = max_val - 2 
            bits_per_entry = int(np.ceil(np.log2(offseted_max_val+1)))
        else:
            bits_per_entry = 0

        data_bitio.write(num_entries, consts.AMAX_NUM_ENTRIES_LEN*8)
        data_bitio.write(bits_per_entry, consts.AMAX_BITS_PER_ENTRY_LEN*8)

        mask_non_ones = self.vector != 1

        for flag, value in zip(mask_non_ones, self.vector):
            data_bitio.write(flag, consts.AMAX_FLAG_BITLEN)

            if flag:
                data_bitio.write(value-2, bits_per_entry)

        return data_bitio
