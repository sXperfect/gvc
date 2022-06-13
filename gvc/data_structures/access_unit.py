from asyncio import constants
from ctypes import util
from dataclasses import dataclass
import logging as log
import typing as t

from gvc.data_structures import block, payload
from . import consts
from .param_set import ParameterSet
from .data_unit import DataUnitHeader
from .block import Block
from ..bitstream import BitstreamReader, RandomAccessHandler
from .. import utils

def comp_blocks_size(blocks):
    size = 0
    for i in range(len(blocks)):
        size += len(blocks[i])
    return size

class AccessUnitHeader(DataUnitHeader):
    access_unit_id:int
    parameter_set_id:int
    num_blocks:int

    def __init__(self, content_len, access_unit_id, parameter_set_id, num_blocks):
        super(AccessUnitHeader, self).__init__(
            consts.DataUnitType.ACCESS_UNIT, 
            content_len
        )

        self.access_unit_id = access_unit_id
        self.parameter_set_id = parameter_set_id
        self.num_blocks = num_blocks

    @classmethod
    def from_blocks(cls, access_unit_id, parameter_set_id, blocks):
        content_len = consts.ACCESS_UNIT_ID_LEN + consts.PARAMETER_SET_ID_LEN + consts.NUM_BLOCKS_LEN
        content_len += comp_blocks_size(blocks)

        num_blocks = len(blocks)

        return cls(
            content_len, access_unit_id, parameter_set_id, num_blocks
        )

    def to_barray(self):

        payload:bytearray = super(AccessUnitHeader, self).to_barray(ret_bytes=False)
        payload += utils.int2bstr(self.access_unit_id, consts.ACCESS_UNIT_ID_LEN)
        payload += utils.int2bstr(self.parameter_set_id, consts.PARAMETER_SET_ID_LEN)
        payload += utils.int2bstr(self.num_blocks, consts.NUM_BLOCKS_LEN)

        return payload

    @classmethod
    def from_bitstream(
        cls, 
        bitstream_reader:BitstreamReader
    ):
        content_len = bitstream_reader.read_bytes(constants.DATA_UNIT_SIZE_LEN)

        access_unit_id = bitstream_reader.read_bytes(consts.ACCESS_UNIT_ID_LEN)
        parameter_set_id = bitstream_reader.read_bytes(consts.PARAMETER_SET_ID_LEN)
        num_blocks = bitstream_reader.read_bytes(consts.NUM_BLOCKS_LEN)

        assert bitstream_reader._byte_aligned(), "Byte not aligned"

        return cls(content_len, access_unit_id, parameter_set_id, num_blocks)


class AccessUnit(object):
    def __init__(
        self,
        header:AccessUnitHeader,
        blocks:t.List
    ):
        # self.access_unit_header = access_unit_header
        self.header = header
        self.blocks = blocks

    @staticmethod
    def blocks_len(blocks):
        size = 0
        for i in range(len(blocks)):
            size += len(blocks[i])

    @property
    def num_blocks(self):
        return len(self.blocks)

    # def get_param_set_id(self):
    #     return self.access_unit_header.parameter_set_id

    # def blocks_to_binary(self):
    #     blocks_b = b''

    #     for block in self.blocks:
    #         blocks_b += block.to_bytes()
        
    #     return blocks_b

    def to_bytes(self):
        payload = self.header.to_barray()
        for i in range(self.num_blocks):
            payload += self.blocks[i].to_bytes()

        return bytes(payload)

    @classmethod
    def from_bitstream(cls, 
        istream:BitstreamReader, 
        parameter_sets:t.List[ParameterSet]
    ):
        start_pos = istream.tell()

        header = AccessUnitHeader.from_bitstream(istream)
        param_set = parameter_sets[header.parameter_set_id]

        blocks = [None] * header.num_blocks
        for i in range(header.num_blocks):
            blocks[i] = Block.from_bitstream(istream, param_set)

        return cls(header, blocks)

    @classmethod
    def from_blocks(
        cls,
        access_unit_id,
        parameter_set_id,
        blocks
    ):
        header = AccessUnitHeader.from_blocks(
            access_unit_id,
            parameter_set_id,
            blocks
        )

        return cls(header, blocks)