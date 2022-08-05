import logging as log
import typing as t
import numpy as np
from ..bitstream import BitIO
import math

def decode_rowcolids(
    data,
    num_entries
):

    ids = np.zeros(num_entries, dtype=np.uint16)
    decode_rowcolids_loop(data, ids, num_entries)

    return ids

def decode_rowcolids_loop(
    data,
    ids,
    num_entries
):

    i_byte = 0
    buffer = 0
    buffer_nbits = 0
    bits_per_id = math.ceil(math.log2(num_entries))

    for i_entry in range(num_entries):        
        while buffer_nbits < bits_per_id:
            buffer <<= 8
            buffer |= data[i_byte]
            
            i_byte += 1
            buffer_nbits += 8
            
        buffer_nbits = buffer_nbits-bits_per_id
        ids[i_entry] = (buffer >> buffer_nbits)
        mask = (1 << buffer_nbits) - 1
        buffer = buffer & mask
        
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
    def from_bytes(cls, data:bytes, num_entries):

        # TODO: Improve by using the cython function
        ids = decode_rowcolids(np.frombuffer(data, dtype=np.uint8), num_entries)
        # data = np.array(np.frombuffer(data, dtype=np.uint8))
        # ids = crc_id.decode_rowcolids(data, num_entries)
        
        return cls(ids)
    
    @classmethod
    def from_randomaccesshandler(cls, ra_handler, num_ids):
        return cls.from_bytes(ra_handler.read(), num_ids)
