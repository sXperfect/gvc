#distutils: language = c++
#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
def cget_col_ids(
    np.ndarray[np.uint32_t, ndim=1, cast=True] qci,
    int p
):

    cdef np.ndarray[np.uint32_t, ndim=1] tqci = np.empty(len(qci)*p, dtype=np.uint32)
    set_tqci(tqci, qci, p, len(qci))

    return tqci

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
cdef void set_tqci(
    np.uint32_t[:] tqci,
    np.uint32_t[:] qci,
    int p,
    int len_qci,
):
    cdef:
        int i,j

    for j in range(len_qci):
        for i in range(p):
            tqci[j*p+i] = qci[j] * p + i