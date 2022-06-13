
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

    def from_bytes(self, ):
        raise NotImplementedError()

    # def header_to_bytes(
    #     self, 
    #     payload_len
    # ):

    #     header_bitio = bitstream.BitIO()
    #     header_bitio.write(consts.DataUnitType.PARAMETER_SET, consts.DATA_UNIT_TYPE_LEN * 8)
    #     header_bitio.write(
    #         consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN + payload_len,
    #         consts.DATA_UNIT_SIZE_LEN*8
    #     )

    #     return header_bitio.to_bytes()
