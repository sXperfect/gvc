import logging as log
import numpy as np
from .utils import bstr2int, int2bstr

class BitstreamWriter(object):
    def __init__(self, f):
        self._accumulator = 0
        self._bit_count = 0
        self.out = f

    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush()

    def __del__(self):
        try:
            self.flush()
        except ValueError:  # I/O operation on closed file
            pass


    def _write_bit(self, bit):
        if self._bit_count == 8:
            self.flush()
        if bit > 0:
            self._accumulator |= 1 << (7 - self._bit_count)
        self._bit_count += 1


    def write_bits(self, bits, nbits):
        if bits.bit_length() > nbits:
            log.error('bit length of bits exceeds nbits')
            raise ValueError('bit length of bits exceeds nbits')

        # Greater than 0 to avoid infinite-loop
        while nbits > 0:
            self._write_bit(bits & 1 << (nbits - 1))
            nbits -= 1


    def byte_aligned(self):
        if self._bit_count % 8 != 0:
            return False
        return True

    def isempty(self):
        return self._bit_count == 0

    def check_empty(self):
        if not self.isempty():
            raise ValueError("bitstream contains data in accumulator")

    def write(self, binary):
        self.check_empty()
        self.out.write(binary)

    def seek(self, offset, whence=0):
        self.check_empty()
        self.out.seek(offset, whence)


    def tell(self):
        self.check_empty()
        return self.out.tell()


    def flush(self):
        # Write only when empty to avoid writing b'\x00'
        if not self.isempty():
            self.out.write(bytearray([self._accumulator]))

        self._accumulator = 0
        self._bit_count = 0


class BitstreamReader(object):
    def __init__(self, f):
        self.input = f
        self._reset()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def _reset(self):
        self._accumulator = 0
        self._bit_count = 0
        self.read = 0

    def _read_bit(self):
        if not self._bit_count:
            a = self.input.read(1)
            
            if a:
                self._accumulator = ord(a)

            self._bit_count = 8
            self.read = len(a)
        rv = (self._accumulator & (1 << (self._bit_count - 1))) >> (self._bit_count - 1)
        self._bit_count -= 1
        return rv

    def _byte_aligned(self):
        if self._bit_count % 8 != 0:
            return False
        return True

    def align_to_byte(self):
        if not self._byte_aligned():
            self.read_bits(self._bit_count)

    def read_bits(self, n):
        v = 0
        for __ in range(n):
            v = (v << 1) | self._read_bit()
        return v
    
    def read_bytes(self, n, ret_int=False):
        assert self._byte_aligned()
        
        payload = self.input.read(n)
        
        if ret_int:
            return bstr2int(payload)
        else:
            return payload

    def seek(self, offset, whence=0):
        self._reset()
        self.input.seek(offset, whence)

    def tell(self):
        if not self._byte_aligned():
            raise ValueError("bitstream must be byte-aligned")
        else:
            if self._bit_count:
                return self.input.tell() - self._bit_count // 8
            else:
                return self.input.tell()

class RandomAccessHandler(object):

    def __init__(self, f:BitstreamReader, start_pos, length, autoseek=False):
        self.f = f
        self.start_pos = start_pos
        self.length = length

        if autoseek:
            self.f.seek(self.length, 1)

    def __call__(self):
        return self.read()
        
    def read(self, nbytes=None):
        self.f.seek(self.start_pos)

        if nbytes is None:
            return self.f.read_bytes(self.length)

        elif nbytes <= self.length:
            return self.f.read_bytes(nbytes)

    # To automate read()
    def __radd__(self, other):
        return other + self.read()

    def __len__(self):
        return self.length

class BitIO(object):
    BYTEORDER = 'big'

    def __init__(self,
        bits=None,
        nbits=None,
    ):

        self.data = 0
        self.len = 0

        # Initialization with initial data
        if bits is not None and nbits is not None:
            self.write(bits, nbits)

        # Initialization without initial data
        elif bits is None and nbits is None:
            pass

        else:
            log.error('Either bits or nbits is None')
            raise RuntimeError("Either bits or nbits is None")

    def write(self, bits, nbits):

        if not isinstance(bits, int):
            bits = int(bits)

        if bits.bit_length() > nbits:
            log.error('bit length of bits exceeds nbits')
            raise RuntimeError("Bit length of bits exceeds nbits")

        self.len += nbits

        self.data <<= nbits
        self.data ^= bits

    def _byte_aligned(self):
        if self.len % 8 == 0:
            return True
        return False

    def len_in_byte(self):
        return int(np.ceil(self.len / 8))

    def __len__(self):
        return self.len

    def align_to_byte(self):
        if not self._byte_aligned():
            self.write(0, 8 - (self.len % 8))

    def to_bytes(self, align=False):

        if align:
            self.align_to_byte()

        if self._byte_aligned():
            return int2bstr(self.data, self.len_in_byte(), order=self.BYTEORDER)
        else:
            raise RuntimeError("Byte not aligned")