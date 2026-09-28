from io import BytesIO

import pytest

from gvc.bitstream import BitIO, BitstreamReader, BitstreamWriter, RandomAccessHandler


def test_bitio_alignment_and_bytes():
    bits = BitIO()
    bits.write(0b101, 3)
    bits.write(0b11, 2)
    assert bits.to_bytes(align=True) == bytes([0b10111000])


def test_writer_and_reader_roundtrip():
    raw = BytesIO()
    writer = BitstreamWriter(raw)
    writer.write_bits(0xA, 4)
    writer.write_bits(0x5, 4)
    writer.flush()

    assert raw.getvalue() == b"\xA5"

    reader = BitstreamReader(BytesIO(raw.getvalue()))
    assert reader.read_bits(4) == 0xA
    assert reader.read_bits(4) == 0x5


def test_read_bytes_requires_alignment():
    reader = BitstreamReader(BytesIO(b"\xF0\x01"))
    assert reader.read_bits(1) == 1
    with pytest.raises(ValueError, match="byte-aligned"):
        reader.read_bytes(1)


def test_random_access_handler_bounds_reads():
    reader = BitstreamReader(BytesIO(b"abcdef"))
    handler = RandomAccessHandler(reader, start_pos=2, length=3)
    assert handler.read() == b"cde"
    assert handler.read(2) == b"cd"
    with pytest.raises(ValueError):
        handler.read(4)



def test_read_bytes_rejects_truncated_input():
    reader = BitstreamReader(BytesIO(b"\x01"))
    with pytest.raises(EOFError, match="unexpected end"):
        reader.read_bytes(2)



def test_read_bits_rejects_end_of_stream():
    reader = BitstreamReader(BytesIO(b""))
    with pytest.raises(EOFError, match="unexpected end"):
        reader.read_bits(1)
