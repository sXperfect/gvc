from . import consts
from .payload import GenotypePayload
from ..bitstream import BitstreamReader
from .. import utils


def _fits_unsigned(value, length_bytes):
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and 0 <= value < (1 << (8 * length_bytes))
    )


class BlockHeader:
    def __init__(self, content_id, block_payload_size):
        if not _fits_unsigned(content_id, consts.CONTENT_ID_LEN):
            raise ValueError("content_id is outside the serialized range")
        if not isinstance(block_payload_size, int) or isinstance(block_payload_size, bool):
            raise TypeError("block payload size must be an integer")
        if block_payload_size < 0:
            raise ValueError("block payload size must be non-negative")
        if not _fits_unsigned(block_payload_size, consts.BLOCK_PAYLOAD_SIZE_LEN):
            raise ValueError("block payload size is outside the serialized range")
        self.content_id = content_id
        self.block_payload_size = block_payload_size

    def to_bytes(self):
        return (
            utils.int2bstr(int(self.content_id), consts.CONTENT_ID_LEN)
            + utils.int2bstr(
                self.block_payload_size, consts.BLOCK_PAYLOAD_SIZE_LEN
            )
        )

    def __len__(self):
        return consts.CONTENT_ID_LEN + consts.BLOCK_PAYLOAD_SIZE_LEN

    @classmethod
    def from_bitstream(cls, bitstream_reader: BitstreamReader):
        content_id = bitstream_reader.read_bytes(
            consts.CONTENT_ID_LEN, ret_int=True
        )
        block_payload_size = bitstream_reader.read_bytes(
            consts.BLOCK_PAYLOAD_SIZE_LEN, ret_int=True
        )
        return cls(content_id, block_payload_size)


class Block:
    def __init__(self, block_header: BlockHeader, block_payload):
        if not isinstance(block_header, BlockHeader):
            raise TypeError("block_header must be a BlockHeader")
        if not isinstance(block_payload, GenotypePayload):
            raise TypeError("block_payload must be a GenotypePayload")
        if block_header.block_payload_size != len(block_payload):
            raise ValueError("block header payload size does not match payload")
        self.block_header = block_header
        self.block_payload = block_payload

    def to_bytes(self):
        return self.block_header.to_bytes() + self.block_payload.to_bytes()

    def __len__(self):
        return len(self.block_header) + len(self.block_payload)

    def stat(self):
        return dict(header=len(self.block_header), **self.block_payload.stat())

    @classmethod
    def from_bitstream(cls, istream, parameter_set):
        block_header = BlockHeader.from_bitstream(istream)
        if not istream._byte_aligned():
            raise ValueError("block header is not byte-aligned")

        if block_header.content_id != consts.ContentID.GENOTYPE:
            raise ValueError(
                "unsupported block content id: {}".format(block_header.content_id)
            )

        block_payload = GenotypePayload.from_bitstream(
            istream,
            parameter_set,
            block_header.block_payload_size,
        )
        return cls(block_header, block_payload)

    @classmethod
    def from_encoded_variant(cls, enc_var: GenotypePayload):
        if not isinstance(enc_var, GenotypePayload):
            raise TypeError("enc_var must be a GenotypePayload")
        block_header = BlockHeader(consts.ContentID.GENOTYPE, len(enc_var))
        return cls(block_header, enc_var)
