import io
import unittest

from gvc.bitstream import BitstreamReader, BitstreamWriter


class TestBitstream(unittest.TestCase):
    def test_bit_roundtrip(self):
        stream = io.BytesIO()
        writer = BitstreamWriter(stream)
        writer.write_bits(0b101, 3)
        writer.write_bits(0b11, 2)
        writer.flush()

        stream.seek(0)
        reader = BitstreamReader(stream)
        self.assertEqual(reader.read_bits(3), 0b101)
        self.assertEqual(reader.read_bits(2), 0b11)

    def test_write_rejects_value_too_wide(self):
        stream = io.BytesIO()
        writer = BitstreamWriter(stream)
        with self.assertRaises(ValueError):
            writer.write_bits(0b100, 2)


if __name__ == "__main__":
    unittest.main()
