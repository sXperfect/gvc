#distutils: language = c++
#cython: language_level=3

import math

import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck.fold(False)
def decode_rowcolids(
    np.ndarray[np.uint8_t, ndim=1, cast=True] data,
    int num_entries
):
    """Decode row/column permutation ids from a read-only byte buffer."""
    if num_entries < 0:
        raise ValueError("num_entries must be non-negative")

    cdef np.ndarray[np.uint16_t, ndim=1] ids = np.zeros(
        num_entries, dtype=np.uint16
    )

    if num_entries <= 1:
        return ids

    cdef int bits_per_id = int(math.ceil(math.log2(num_entries)))
    cdef int required_bytes = (num_entries * bits_per_id + 7) // 8
    if data.shape[0] < required_bytes:
        raise ValueError("permutation payload is truncated")

    decode_rowcolids_loop(data, ids, num_entries, bits_per_id)

    if np.any(ids >= num_entries):
        raise ValueError("permutation id is outside the valid range")

    return ids


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck.fold(False)
cdef void decode_rowcolids_loop(
    const np.uint8_t[:] data,
    np.uint16_t[:] ids,
    int num_entries,
    int bits_per_id
):
    cdef:
        int i_entry
        int delta_nbits
        int i_byte = 0
        unsigned int buffer = 0
        int buffer_nbits = 0
        unsigned int mask

    for i_entry in range(num_entries):
        while buffer_nbits < bits_per_id:
            buffer <<= 8
            buffer |= data[i_byte]
            i_byte += 1
            buffer_nbits += 8

        delta_nbits = buffer_nbits - bits_per_id
        buffer_nbits = delta_nbits
        ids[i_entry] = buffer >> delta_nbits
        mask = (1 << delta_nbits) - 1
        buffer &= mask
