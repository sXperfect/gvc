import typing as t

from . import consts
from .block import Block
from .data_unit import DataUnitHeader
from ..bitstream import BitstreamReader
from .. import utils


def _fits_unsigned(value, length_bytes):
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and 0 <= value < (1 << (8 * length_bytes))
    )


def comp_blocks_size(blocks):
    return sum(len(block) for block in blocks)


class AccessUnitHeader(DataUnitHeader):
    def __init__(self, content_len, access_unit_id, parameter_set_id, num_blocks):
        fixed_content_len = (
            consts.ACCESS_UNIT_ID_LEN
            + consts.PARAMETER_SET_ID_LEN
            + consts.NUM_BLOCKS_LEN
        )
        if content_len < fixed_content_len:
            raise ValueError("access-unit content length is smaller than its header")
        if not _fits_unsigned(access_unit_id, consts.ACCESS_UNIT_ID_LEN):
            raise ValueError("access_unit_id is outside the serialized range")
        if not _fits_unsigned(parameter_set_id, consts.PARAMETER_SET_ID_LEN):
            raise ValueError("parameter_set_id is outside the serialized range")
        if not _fits_unsigned(num_blocks, consts.NUM_BLOCKS_LEN):
            raise ValueError("num_blocks is outside the serialized range")

        super().__init__(consts.DataUnitType.ACCESS_UNIT, content_len)
        self.access_unit_id = access_unit_id
        self.parameter_set_id = parameter_set_id
        self.num_blocks = num_blocks

    @classmethod
    def from_blocks(cls, access_unit_id, parameter_set_id, blocks):
        blocks = list(blocks)
        content_len = (
            consts.ACCESS_UNIT_ID_LEN
            + consts.PARAMETER_SET_ID_LEN
            + consts.NUM_BLOCKS_LEN
            + comp_blocks_size(blocks)
        )
        return cls(content_len, access_unit_id, parameter_set_id, len(blocks))

    def to_barray(self):
        payload = super().to_barray(ret_bytes=False)
        payload += utils.int2bstr(self.access_unit_id, consts.ACCESS_UNIT_ID_LEN)
        payload += utils.int2bstr(self.parameter_set_id, consts.PARAMETER_SET_ID_LEN)
        payload += utils.int2bstr(self.num_blocks, consts.NUM_BLOCKS_LEN)
        return payload

    @classmethod
    def from_bitstream(cls, bitstream_reader: BitstreamReader):
        total_len = bitstream_reader.read_bytes(consts.DATA_UNIT_SIZE_LEN, ret_int=True)
        base_header_len = consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN
        fixed_content_len = (
            consts.ACCESS_UNIT_ID_LEN
            + consts.PARAMETER_SET_ID_LEN
            + consts.NUM_BLOCKS_LEN
        )
        if total_len < base_header_len + fixed_content_len:
            raise ValueError("access-unit length is smaller than its header")

        content_len = total_len - base_header_len
        access_unit_id = bitstream_reader.read_bytes(
            consts.ACCESS_UNIT_ID_LEN, ret_int=True
        )
        parameter_set_id = bitstream_reader.read_bytes(
            consts.PARAMETER_SET_ID_LEN, ret_int=True
        )
        num_blocks = bitstream_reader.read_bytes(consts.NUM_BLOCKS_LEN, ret_int=True)

        if not bitstream_reader._byte_aligned():
            raise ValueError("access-unit header is not byte-aligned")

        return cls(content_len, access_unit_id, parameter_set_id, num_blocks)


class AccessUnit:
    def __init__(self, header: AccessUnitHeader, blocks: t.List[Block]):
        if not isinstance(header, AccessUnitHeader):
            raise TypeError("header must be an AccessUnitHeader")
        blocks = list(blocks)
        if any(not isinstance(block, Block) for block in blocks):
            raise TypeError("access-unit blocks must be Block instances")
        if header.num_blocks != len(blocks):
            raise ValueError("access-unit header block count does not match payload")
        self.header = header
        self.blocks = blocks

    @staticmethod
    def blocks_len(blocks):
        return comp_blocks_size(blocks)

    @property
    def num_blocks(self):
        return len(self.blocks)

    def to_bytes(self):
        payload = self.header.to_barray()
        for block in self.blocks:
            payload += block.to_bytes()
        return bytes(payload)

    @classmethod
    def from_bitstream(cls, istream: BitstreamReader, parameter_sets):
        start_pos = istream.tell()
        header = AccessUnitHeader.from_bitstream(istream)

        try:
            param_set = parameter_sets[header.parameter_set_id]
        except (KeyError, IndexError) as exc:
            raise ValueError(
                "missing parameter set {}".format(header.parameter_set_id)
            ) from exc

        blocks = [Block.from_bitstream(istream, param_set) for _ in range(header.num_blocks)]

        expected_consumed = consts.DATA_UNIT_SIZE_LEN + header.content_len
        consumed = istream.tell() - start_pos
        if consumed != expected_consumed:
            raise ValueError(
                "access-unit length mismatch: expected {}, consumed {}".format(
                    expected_consumed, consumed
                )
            )

        return cls(header, blocks)

    @classmethod
    def from_blocks(cls, access_unit_id, parameter_set_id, blocks):
        blocks = list(blocks)
        header = AccessUnitHeader.from_blocks(access_unit_id, parameter_set_id, blocks)
        return cls(header, blocks)
