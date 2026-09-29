from dataclasses import dataclass

from . import consts
from .. import utils


def _fits_unsigned(value, length_bytes):
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and 0 <= value < (1 << (8 * length_bytes))
    )


@dataclass
class DataUnitHeader:
    type: int
    content_len: int

    def __post_init__(self):
        if not _fits_unsigned(self.type, consts.DATA_UNIT_TYPE_LEN):
            raise ValueError("data-unit type is outside the serialized range")
        if not isinstance(self.content_len, int) or isinstance(self.content_len, bool):
            raise TypeError("data-unit content length must be an integer")
        if self.content_len < 0:
            raise ValueError("data-unit content length must be non-negative")

        max_total = (1 << (8 * consts.DATA_UNIT_SIZE_LEN)) - 1
        if len(self) > max_total:
            raise ValueError("data-unit length is outside the serialized range")

    def to_barray(self, ret_bytes=False):
        payload = bytearray()
        payload += utils.int2bstr(int(self.type), consts.DATA_UNIT_TYPE_LEN)
        payload += utils.int2bstr(len(self), consts.DATA_UNIT_SIZE_LEN)
        return bytes(payload) if ret_bytes else payload

    def __len__(self):
        return consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN + self.content_len

    @classmethod
    def from_bitstream(cls, data_unit_type, bitreader):
        total_len = bitreader.read_bytes(consts.DATA_UNIT_SIZE_LEN, ret_int=True)
        header_len = consts.DATA_UNIT_TYPE_LEN + consts.DATA_UNIT_SIZE_LEN
        if total_len < header_len:
            raise ValueError("data-unit length is smaller than its header")
        return cls(data_unit_type, total_len - header_len)
