
from dataclasses import dataclass
import abc
import logging as log
import typing as t
from . import consts
from .. import bitstream
from .. import utils

@dataclass
class DataUnitHeader():
    type:int
    content_len:int # Length of the content

    def to_barray(self, ret_bytes=False):
        payload = bytearray()
        payload += utils.int2bstr(self.type, consts.DATA_UNIT_TYPE_LEN)
        payload += utils.int2bstr(
            len(self),
            consts.DATA_UNIT_SIZE_LEN
        )

        return payload

    def __len__(self):
        return consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN + self.content_len

    # def from_bytes(self, ):
    #     raise NotImplementedError()
    
    @classmethod
    def from_bitstream(cls, type, bitreader):
        
        total_len = bitreader.read_bytes(consts.DATA_UNIT_SIZE_LEN, ret_int=True)
        content_len = total_len - (consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN)
        
        return cls(type, content_len)