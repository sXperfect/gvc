#distutils: language = c++
#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython
import math

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
def decode_rowcolids(
    np.ndarray[np.uint8_t, ndim=1, cast=True] data,
    int num_entries
):

    cdef np.ndarray[np.uint16_t, ndim=1] ids = np.zeros(num_entries, dtype=np.uint16)
    decode_rowcolids_loop(data, ids, num_entries)

    return ids

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
cdef void decode_rowcolids_loop(
    np.uint8_t[:] data,
    np.uint16_t[:] ids,
    int num_entries
):

    cdef:
        int i_entry
        np.int8_t delta_nbits
        np.uint16_t i_byte = 0
        np.uint16_t buffer = 0
        np.int8_t buffer_nbits = 0
        np.uint8_t bits_per_id = math.ceil(math.log2(num_entries))

    for i_entry in range(num_entries):        
        while buffer_nbits < bits_per_id:
            buffer <<= 8
            buffer |= data[i_byte]
            
            i_byte += 1
            buffer_nbits += 8
            
        delta_nbits = buffer_nbits-bits_per_id
        buffer_nbits = delta_nbits
        ids[i_entry] = (buffer >> delta_nbits)
        mask = (1 << delta_nbits) - 1
        buffer = buffer & mask
        