#distutils: language = c++
#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free
from libc.string cimport memset
from libcpp.string cimport string

ctypedef np.npy_bool bool
ctypedef np.uint8_t uint8
ctypedef np.int8_t int8
# ctypedef object* obj_ptr

cdef extern from "Python.h":
    object PyBytes_FromStringAndSize(char *s, Py_ssize_t length)
    # object PyString_FromStringAndSize(char *s, Py_ssize_t length)
    object PyLong_FromSsize_t(Py_ssize_t v)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
def debin_rc_bin_split(np.ndarray[bool, ndim=2, cast=True] bin_mat, np.ndarray[uint8, ndim=1] bitlen_vect):

    cdef int nrows = bitlen_vect.shape[0]
    cdef int ncols = bin_mat.shape[1]

    cdef np.ndarray[np.uint8_t, ndim=2] mat = np.zeros((nrows, ncols), dtype=np.uint8)

    debin_rc_bin_split_loop(mat, bin_mat, bitlen_vect, nrows, ncols)

    return mat

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
cdef void debin_rc_bin_split_loop(
    np.uint8_t[:, :] mat,
    # np.uint8_t[:, :] bin_mat,
    bool[:, :] bin_mat,
    np.uint8_t[:] bitlen_vect,
    int nrows,
    int ncols
):

    cdef:
        int i, j, k
        int irow_bin_mat = 0

    for i in range(nrows):
        for j in range(bitlen_vect[i]):
            for k in range(ncols):
                mat[i, k] |= (<uint8>bin_mat[irow_bin_mat+j, k]) << j

        irow_bin_mat += bitlen_vect[i]

cdef char gt_val_to_gt_char(np.int8_t v):

    if v == -1:
        return 46 # '.'
    else:
        return v+48 # 0, 1, ..., 9

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
def recon_gt_mat_with_phase_val(np.ndarray[np.int8_t, ndim=2] allele_mat, bool phase_val, int p):

    cdef Py_ssize_t length
    cdef char* gt_mat

    crecon_gt_mat_with_phase_val(
        &gt_mat, &length, allele_mat, phase_val, p
    )

    out_gt_mat = PyBytes_FromStringAndSize(gt_mat, length)
    free(gt_mat)

    return out_gt_mat

@cython.boundscheck(False) #? turn off bounds-checking for entire function
@cython.wraparound(False)  #? turn off negative index wrapping for entire function
@cython.overflowcheck.fold(False)
cdef void crecon_gt_mat_with_phase_val(
    char** gt_mat, 
    Py_ssize_t *length, 
    # np.ndarray[np.uint8_t, ndim=2] allele_mat, 
    np.int8_t[:, :] allele_mat,
    bool phase_val, 
    int p
):

    cdef:
        int n_cols = allele_mat.shape[1]
        int n_variants = allele_mat.shape[0]
        int n_samples = n_cols//p
        char phase_char = b'/'

        # int max_val = np.amax(allele_mat)
        int char_per_genotype = 0

        Py_ssize_t gt_mat_len = 0

        Py_ssize_t n = 0
        int8 curr_gt_val = 0
        int8 variant_max_val = 0

        int i_sample, k

    char_per_genotype = 1

    #? length: number of genotypes{., 0, 1, ...} + number of phasing{\, |}  + number of separator {\t, \n}
    #? Multiplied by maximum genotype value (e.g 10 requires 2 bytes)
    gt_mat_len = (n_cols*char_per_genotype + (p-1)*n_samples + n_samples) * n_variants

    if phase_val == 0:
        phase_char = b'|'

    gt_mat[0] = <char *> malloc((gt_mat_len + 1) * sizeof(char))
        
    #? Needed only if "PyBytes_FromStringAndSize" is not used in the next step
    #? This is used to avoid string remnant from other process
    #? memset(gt_mat[0], b'\0', (gt_mat_len + 1) * sizeof(char))
    memset(gt_mat[0], b'0', (gt_mat_len + 1) * sizeof(char))

    for i_variant in range(n_variants):
        for i_sample in range(n_samples):

            curr_gt_val = allele_mat[i_variant, i_sample*p]
            
            #? Handle case where GT for current sample is available
            if curr_gt_val == -2:
                continue

            if curr_gt_val != 0:
                gt_mat[0][n] = gt_val_to_gt_char(curr_gt_val)
            n += 1

            for k in range(1, p):
                curr_gt_val = allele_mat[i_variant, i_sample*p+k]
                if curr_gt_val == -2:
                    break

                gt_mat[0][n] = phase_char
                if curr_gt_val != 0:
                    gt_mat[0][n] = gt_val_to_gt_char(curr_gt_val)
                n += 2
            
            #? Add separator at the end of processing one sample
            if i_sample < n_samples-1:
                gt_mat[0][n] = b'\t'
            else:
                gt_mat[0][n] = b'\n'
            n += 1

    length[0] = n