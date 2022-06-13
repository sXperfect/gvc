import logging as log
import typing as t
import numpy as np
from . import consts
from ..bitstream import BitIO
from .. import libds

class RowColIds(object):
    def __init__(self, ids):
        self.ids = ids

    def to_bitio(self) -> BitIO:
        num_entries = self.ids.shape[0]
        bits_per_id = np.ceil(np.log2(num_entries)).astype(int)

        ids_bitio = BitIO()

        for curr_id in self.ids:
            ids_bitio.write(int(curr_id), int(bits_per_id))

        return ids_bitio

    def to_bytes(self):
        return self.to_bitio().to_bytes(align=True)

    @classmethod
    def from_bytes(cls, data:bytes, num_ids):
        # data_bytes = io.BytesIO(data)
        # istream = gvc.bitstream.BitstreamReader(data_bytes)

        # bits_per_id = np.ceil(np.log2(num_ids)).astype(int)

        # ids = np.zeros((num_ids),dtype=int)
        # for i_entry in range(num_ids):
        #     ids[i_entry] = istream.read_bits(bits_per_id)

        # return cls(ids)

        return cls(
            libds.decode_rowcolids(
                data, num_ids
            )
        )

    @classmethod
    def from_randomaccesshandler(cls, ra_handler, num_ids):
        return cls.from_bytes(ra_handler.read(), num_ids)
