from typing import List
import logging as log
import time

import numpy as np
import itertools as it

import gvc.common

from .data_structures.consts import BinarizationID
from . import cdebinarize
from . import debinarize

_phasing_dict = {
    '|' : 0,
    '/' : 1
}

_phasing_val_to_str = ['|', '/']

def _gt_code_to_int(gt_code):
    try:
        return gvc.common.ALLELE_DTYPE(gt_code)
    except ValueError:
        if gt_code == ".":
            return gvc.common.SIGNED_ALLELE_DTYPE(-1)
        else:
            raise ValueError()

def _int_to_gt_code(val):

    if val >= 0:
        return str(val)
    else:
        return "."

def matrix_to_tensor(matrix, num_matrix):

    list_matrix = np.split(
        np.expand_dims(matrix, axis=1), 
        matrix.shape[1]//num_matrix, 
        axis=2
    )

    return np.concatenate(list_matrix, axis=1)

def _tensor_to_matrix(tensor):
    """
    Transform tensor to matrix.
    Part of the implementation 4.1

    Parameters
    ----------
    tensor : ndarray
        a ndarray with dimension of 3

    Returns
    -------
    matrix : ndarray
        a ndarray with dimension of 2
    """
    assert(len(tensor.shape) == 3)

    list_matrix = np.split(tensor, tensor.shape[1], axis=1)

    matrix = np.concatenate(list_matrix, axis=2).squeeze(axis=1)
    
    assert(np.issubdtype(matrix.dtype, tensor.dtype))
    
    return matrix

def store_tensors_to_file(f, allele_tensor, phasing_tensor):
    m, n, p = allele_tensor.shape

    for i in range(m):

        txt = [None] * n

        for j in range(n):

            sample_val = ""

            if allele_tensor[i, j, 0] != -2:
                allele_val = allele_tensor[i, j, 0]
                sample_val += _int_to_gt_code(allele_val)

                for k in range(1, p):
                    allele_val = allele_tensor[i, j, k]
                    if allele_val == -2:
                        break

                    phasing_val = phasing_tensor[i, j, k-1]

                    sample_val += _phasing_val_to_str[phasing_val]
                    sample_val += _int_to_gt_code(allele_val)

                txt[j] = sample_val
            else:
                txt[j] = ""

        f.write('\t'.join(txt) + '\n')


def tensor_to_txt(allele_tensor, phasing_tensor):
    m, n, p = allele_tensor.shape

    max_val = allele_tensor.max()
    gt_phased  =  np.array(["|".join(x) for x in list(it.product(np.arange(max_val+1).astype(str), repeat=p))])
    gt_unphased = np.array(["/".join(x) for x in list(it.product(np.arange(max_val+1).astype(str), repeat=p))])

    # lines = np.tile(gt_phased[0], (m,n))
    lines = np.tile(gt_phased[0], (m,n))

    mask = (allele_tensor >= 0).all(axis=2)

    phased_mask = (phasing_tensor == 0).all(axis=2) & mask
    phased_idx = allele_tensor[phased_mask].dot((max_val+1)**np.flip(np.arange(p)))
    lines[phased_mask] = gt_phased[phased_idx]

    unphased_mask = (phasing_tensor == 1).all(axis=2) & mask
    unphased_idx = allele_tensor[unphased_mask].dot((max_val+1)**np.flip(np.arange(p)))
    lines[unphased_mask] = gt_unphased[unphased_idx]

    rest_mask = ~(phased_mask | unphased_mask)
    
    if rest_mask.any():

        for i in range(m):
            for j in range(n):
                if rest_mask[i,j]:
                    sample_val = ""
                    if allele_tensor[i, j, 0] != -2:
                        allele_val = allele_tensor[i, j, 0]
                        sample_val += _int_to_gt_code(allele_val)

                        for k in range(1, p):
                            allele_val = allele_tensor[i, j, k]
                            if allele_val == -2:
                                break

                            phasing_val = phasing_tensor[i, j, k-1]

                            sample_val += _phasing_val_to_str[phasing_val]
                            sample_val += _int_to_gt_code(allele_val)

                        lines[i,j] = sample_val
                    else:
                        lines[i,j] = ""

    return lines

# def simd_tensor_to_txt(allele_tensor, phasing_tensor):
#     m, n, p = allele_tensor.shape

#     max_val = allele_tensor.max()
#     avail_allele_vals = '.' + "".join(np.arange(max_val+1).astype(str))
#     avail_phase_val = '|/'

#     if p == 1:
#         codebook = avail_allele_vals
#     elif p > 1:
#         codebook = ["".join(l) for l in it.product(avail_allele_vals, avail_phase_val, repeat=p-1)]
#         codebook = ["".join(l) for l in it.product(codebook, avail_allele_vals)]

#     def_val_idx = 0
#     for k in range(1, p):
#         def_val_idx += len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k)
#         # def_val_idx += 0*len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k-1)
#     def_val_idx += 1

#     allele_tensor += 1

#     flag_mask = allele_tensor == -1
#     complete_allele_mask = (np.logical_not(flag_mask)).all(axis=2)

#     # gt_matrix_str = np.full((m,n), codebook[def_val_idx])
#     gt_matrix_str = np.empty((m,n), dtype='<U{}'.format(2*p))

#     index_mat = np.zeros(gt_matrix_str.shape, dtype=np.int)
#     for k in range(1, p):
#         allele_factor = len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k)
#         index_mat[complete_allele_mask] += allele_tensor[complete_allele_mask, k-1] * allele_factor
#         phase_factor = len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k-1)
#         index_mat[complete_allele_mask] += phasing_tensor[complete_allele_mask, k-1] * phase_factor

#     index_mat[complete_allele_mask] += allele_tensor[complete_allele_mask, -1]
#     gt_matrix_str[complete_allele_mask] = np.array(codebook)[index_mat[complete_allele_mask]]

#     if (flag_mask).any():
#         for i, j, k in np.argwhere(flag_mask):
#             sample_val = avail_allele_vals[allele_tensor[i, j, 0]]
#             for l in range(1, k):
#                 sample_val = avail_allele_vals[allele_tensor[i, j, l]] + avail_phase_val[phasing_tensor[i, j, l-1]]

#             gt_matrix_str[i, j] = sample_val

#     delim_mat_str = np.full((m,n,1), '\t')
#     delim_mat_str[:, -1] = '\n'
#     gt_matrix_str = np.concatenate((np.expand_dims(gt_matrix_str,2), delim_mat_str), axis=2)
#     gt_matrix_str = _tensor_to_matrix(gt_matrix_str)
    
#     return "".join(gt_matrix_str.reshape(-1))

def compare_tensor_to_txt(allele_tensor, phasing_tensor, f):
    m, n, p = allele_tensor.shape
    start_time = time.time()

    for i in range(m):

        txt = [None] * n

        for j in range(n):

            sample_val = ""

            if allele_tensor[i, j, 0] != -2:
                allele_val = allele_tensor[i, j, 0]
                sample_val += _int_to_gt_code(allele_val)

                for k in range(1, p):
                    allele_val = allele_tensor[i, j, k]
                    if allele_val == -2:
                        break

                    phasing_val = phasing_tensor[i, j, k-1]

                    sample_val += _phasing_val_to_str[phasing_val]
                    sample_val += _int_to_gt_code(allele_val)

                txt[j] = sample_val
            else:
                txt[j] = ""

        recon_line = '\t'.join(txt) + '\n'

        if i % 1024 == 0:
            end_time = time.time()
            log.info("iline {} - {:.02f}s".format(i, end_time - start_time))
            start_time = end_time

        orig_line = f.readline()
        if recon_line != orig_line:
            raise ValueError('Line {}'.format(i))


def reconstruct_genotype_matrix(allele_matrix, phasing_matrix, p):

    allele_tensor = matrix_to_tensor(allele_matrix, p)

    if p >= 1:
        if isinstance(phasing_matrix, int):
            allele_tensor_shape = allele_tensor.shape

            phasing_tensor = np.zeros((*allele_tensor_shape[:-1], p-1), dtype=np.bool)
            phasing_tensor[:, :, :] = phasing_matrix

        else:
            phasing_tensor = matrix_to_tensor(phasing_matrix, p-1)

    else:
        log.error('Invalid value p: {}'.format(p))
        raise ValueError('Invalid value p: {}'.format(p))

    return allele_tensor, phasing_tensor

def adaptive_max_value(allele_matrix):

    allele_matrix = allele_matrix.copy()

    assert(np.issubdtype(allele_matrix.dtype, np.signedinteger))

    dot_mask = allele_matrix == -1
    na_mask = allele_matrix == -2

    if np.any(dot_mask):
        missing_rep_val = np.max(allele_matrix) + 1
        allele_matrix[dot_mask] = missing_rep_val
    else:
        missing_rep_val = None

    if np.any(na_mask):
        na_rep_val = np.max(allele_matrix) + 1
        allele_matrix[na_mask] = na_rep_val
    else:
        na_rep_val = None

    return allele_matrix.astype(gvc.common.ALLELE_DTYPE), missing_rep_val, na_rep_val

def undo_adaptive_max_value(allele_array, missing_rep_val, na_rep_val):
    allele_array = allele_array.astype(gvc.common.SIGNED_ALLELE_DTYPE)

    if na_rep_val is not None:
        allele_array[allele_array == na_rep_val] = -2

    if missing_rep_val is not None:
        allele_array[allele_array == missing_rep_val] = -1
    
    return allele_array

def split_genotype_matrix(genotype_matrix: List[str]):
    r"""
    Split genotype matrix into the corresponding allele and phasing matrices.
    Implements 4.1

    Parameters
    ----------
    genotype_matrix : list of string
        The genotype matrix

    Returns
    -------
    allele_tensor : ndarray
        The allele tensor (3d)
    phasing_tensor : ndarray
        The phasing tensor (3d)
    """

    log.debug('splitting genotype matrix')

    # Convert to numpy array and split by tabulator
    genotype_matrix = np.array([x.split('\t') for x in genotype_matrix])
    genotype_matrix[:,-1] = [x.strip() for x in genotype_matrix[:,-1]] # remove "\n" in last column

    m,n = genotype_matrix.shape
    #p = int((len(genotype_matrix[0][0].split('\t')[0]) + 1) / 2)  # 1 for phasing & 1 for genotype
    p_candidates = ((np.vectorize(len)(genotype_matrix[0,:]) + 1) / 2).astype(int) # TODO: is it sufficient to check only first column?
    p = np.max(p_candidates)

    allele_tensor = np.zeros([m, n, p], dtype=gvc.common.SIGNED_ALLELE_DTYPE)
    assert(np.issubdtype(allele_tensor.dtype, np.signedinteger))

    if p-1 > 0:
        phasing_tensor = np.zeros([m, n, (p - 1)], dtype=gvc.common.PHASING_DTYPE)
    else:
        phasing_tensor = None

    # allele_tensor and phasing_tensor are initialized assuming the default case "0|...0" in genotype_matrix
    all_zeros_template = "|".join(np.repeat('0', p)) # default template

    # Here we construct templates with the most common (binary) cases
    # either all phased or all unphased (we handle the rest later in a slow for loop)
    na_templates_str = np.array(['0', '1'])
    binary_phased   = np.array(["|".join(x) for x in list(it.product('01', repeat=p))])
    binary_phased   = np.concatenate((binary_phased, na_templates_str)) # add single values as well
    binary_unphased = np.array(["/".join(x) for x in list(it.product('01', repeat=p))])
    binary_templates_str = np.concatenate((binary_phased, binary_unphased))

    # Same templates, now as int
    na_templates_int = (-2)*np.ones((2, p))
    na_templates_int[0,0] = 0; na_templates_int[1,0] = 1
    binary_templates_int = np.arange(2**(p))[:, np.newaxis] >> np.arange(p)[::-1] & 1
    binary_templates_int = np.concatenate((binary_templates_int, na_templates_int, binary_templates_int))
    
    # Store cases which we have to handle later
    handle_later_mask = genotype_matrix != all_zeros_template

    # Loop over templates
    for idx, bt in enumerate(binary_templates_str[1:],1):
        mask = genotype_matrix == bt
        allele_tensor[mask] = binary_templates_int[idx,:]
        phasing_tensor[mask] = idx >= len(binary_phased)
        handle_later_mask[mask] = False # We don't need to handle these cases later!


    # Slow loop for more complicated cases only
    if handle_later_mask.any():
        for i in range(m):
            for j in range(n):
                if not handle_later_mask[i,j]:
                    continue # skip # TODO: this is still slow, reduce for loop to iterate only over elements for which mask[i,j] is True
                else:
                    try:
                        allele_tensor[i, j, 0] = _gt_code_to_int(genotype_matrix[i,j][0])
                    except ValueError:
                        log.error('Could not parse row {}: {}'.format(i, genotype_matrix[i,:]))
                        raise ValueError('Could not parse row {}: {}'.format(i, genotype_matrix[i,:]))

                    p_new = (len(genotype_matrix[i,j])+1)//2
                    if p_new > p:
                        # Insert additional depth with values -2, which represent *NotAvailable*
                        allele_tensor = np.concatenate((
                            allele_tensor,
                            -2*np.ones([m, n, p_new - p], dtype=gvc.common.SIGNED_ALLELE_DTYPE)
                        ), axis=2)

                        assert(np.issubdtype(allele_tensor.dtype, np.signedinteger))

                        additional_depth = np.zeros([m, n, p_new-p], dtype=bool)
                        if phasing_tensor is None:
                            phasing_tensor = additional_depth
                        else:
                            phasing_tensor = np.concatenate((phasing_tensor, additional_depth), axis=2)

                        p = p_new

                    for k in range(1, p):
                        try:
                            gt_code = genotype_matrix[i,j][2*k]
                            phase = genotype_matrix[i,j][2*k-1]

                            phasing_tensor[i, j, k-1] = _phasing_dict[phase]
                        except IndexError:
                            ### Adaptive Max Value preprocessing
                            # GT code is set to -2
                            gt_code = -2
                            # Phase of *NotAvailable* is arbitrarily set
                            phase = 0

                            phasing_tensor[i, j, k-1] = phase

                        try:
                            if gt_code != -2:
                                allele_tensor[i, j, k] = _gt_code_to_int(gt_code)
                            else:
                                allele_tensor[i, j, k] = -2
                        except ValueError:
                            log.error('Could not parse row {}: {}'.format(i, genotype_matrix[i,:]))
                            raise ValueError('Could not parse row {}: {}'.format(i, genotype_matrix[i,:]))

    assert(np.issubdtype(allele_tensor.dtype, np.signedinteger))

    if phasing_tensor is not None:
        return _tensor_to_matrix(allele_tensor), _tensor_to_matrix(phasing_tensor), p
    else:
        return _tensor_to_matrix(allele_tensor), None, p


def bin_bit_plane(matrix, axis=None, **kwargs):
    """
    Binarize a matrix into bitplanes.
    Implements 4.2.1

    Parameters
    ----------
    matrix : ndarray
        The matrix to binarize. Must contain integers greater or equal to zero.
    axis : int
        The direction, which matrices are concatenated

    Returns
    -------
    bin_matrix : list of ndarray
        Binarized matrix using bit plane method
    """
    log.debug('binarizing allele matrix using bit plane')

    bit_planes = []

    bit_depth = np.ceil(np.log2(matrix.max() + 1)).astype(gvc.common.ALLELE_DTYPE)

    for i_bit in range(bit_depth):
        bit_tensor = np.bitwise_and(matrix, int(2**i_bit)).astype(gvc.common.BIN_DTYPE)
        bit_planes.append(bit_tensor)
    
    # Concatenate binary matrices
    if len(bit_planes) > 1:
        if axis in (0, 1):
            bin_matrices = [np.concatenate(bit_planes, axis=axis)]
        elif axis == 2:
            bin_matrices = bit_planes
        else:
            log.error('Invalid axis: {}'.format(axis))
            raise gvc.errors.GvcError()
    else:
        bin_matrices = bit_planes

    return bin_matrices, bit_depth

def debin_bit_plane(bin_matrices, bit_depth, axis):

    if axis < 2:
        assert len(bin_matrices) == 1

        bit_planes = np.split(bin_matrices[0], bit_depth, axis)

    elif axis == 2 :
        bit_planes = bin_matrices

    else:
        log.error('Invalid axis: {}'.format(axis))
        raise gvc.errors.GvcError()

    assert bit_depth == len(bit_planes)

    matrix = bit_planes[0].astype(gvc.common.ALLELE_DTYPE)
    for i in range(1, bit_depth):
        matrix += 2**i * bit_planes[i].astype(gvc.common.ALLELE_DTYPE)

    return matrix

def bin_row_bin_split(matrix, **kwargs):
    """
    Binarize a matrix by row splitting.
    Implements 4.2.2

    Parameters
    ----------
    matrix : ndarray
        The matrix to binarize. Must contain integers greater or equal to zero.

    Returns
    -------
    bin_mat : ndarray (2d)
        Binarized matrix
    amax_vect : ndarray (1d)
        Maximum value on each row 
    """

    log.debug('binarizing allele matrix by row splitting')

    nrow = matrix.shape[0]

    bitlen_vect = np.max(matrix, axis=1).astype(np.uint)
    # Force max value 0 to 1
    bitlen_vect[bitlen_vect == 0] = 1

    # consider the case when a line has all 0s
    log.debug("Greatest value in the block: {}".format(bitlen_vect.max()))

    bitlen_vect = np.ceil(np.log2(bitlen_vect + 1)).astype(np.uint16)
    bin_mat_nrows = np.sum(bitlen_vect).astype(np.uint)
    bin_mat_shape = (bin_mat_nrows, matrix.shape[1])
    bin_mat = np.zeros(bin_mat_shape, dtype=np.bool)

    i_b = 0  # row id in bin_mat
    for i in range(nrow):  # iterate over rows in matrix

        bitlen = bitlen_vect[i]

        for i_bit in range(bitlen):
            bin_mat[i_b + i_bit] = np.bitwise_and(matrix[i, :], int(2**i_bit)).astype(np.bool)

        i_b += bitlen

    return bin_mat, bitlen_vect

def debin_row_bin_split(bin_matrices, bitlen_vect, **kwargs):

    try:
        if isinstance(bin_matrices, List):
            bin_mat = bin_matrices[0]
        elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 1:
            bin_mat = bin_matrices[0]
        # elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 3:
        #     bin_mat = bin_matrices[0, :, :]
        # elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 2:
        #     bin_mat = bin_matrices
        else:
            raise TypeError
        
    except TypeError:
        error_msg = "Invalid data type:{}".format(type(bin_matrices))
        log.info(error_msg)
        raise ValueError(error_msg)

    # TODO: fix cython function debin_rc_bin_split
    # matrix = cdebinarize.debin_rc_bin_split(bin_mat, bitlen_vect.astype(np.uint8))
    matrix = debinarize.debin_rc_bin_split(bin_mat, bitlen_vect.astype(np.uint8))

    return matrix

BINARIZATIONS = {
    BinarizationID.BIT_PLANE: {
        "name": "bit_plane",
        "transform": bin_bit_plane,
        "inverse_transform": debin_bit_plane,
    },
    BinarizationID.ROW_BIN_SPLIT: {
        "name": "row_bin_split",
        "transform": bin_row_bin_split,
        "inverse_transform": debin_row_bin_split,
    }
}

AVAIL_BINARIZATION_MODES = [v['name'] for v in BINARIZATIONS.values()]
BINARIZATION_STR2ID = {v["name"]:k for k, v in BINARIZATIONS.items()}

def debinarize_bin_matrices(
    bin_matrices,
    additional_info,
    axis,
    binarization_mode:int
):
    try:
        debinarizer = BINARIZATIONS[binarization_mode]["inverse_transform"]
    except IndexError:
        log.error('Invalid binarization method')
        raise KeyError('Invalid binarization method')

    return debinarizer(bin_matrices, additional_info, axis=axis)

def binarize_allele_matrix(
    matrix, 
    binarization_id:int,
    axis:int=None
):
    """
    Binarize alelle matrix
    Implements 4.2

    Parameters
    ----------
    matrix : ndarray
        The matrix to binarize. Must contains integers greater or equal to zero.
    binarizer : string, optional
        Either 'bit_plane' for binarization using bit plane or 'row_split' by using row splitting method.
        Mode with value None is not allowed
    axis : int, optional
        Direction, in which matrices are concatenated when using bit plane method.
        0 for row direction and 1 for column direction

    Returns:
    ----------
    bin_matrices : list of ndarray
        Binarized allele matrices.
        Number of matrices depends on binarization method.
    additional_info : int or ndarray
        Vector which element describe maximum value of each row of original allele matrix.
        The value is None if mode is not 'row_split'
    """

    log.debug('binarizing allele matrix')
    
    try:
        binarizer = BINARIZATIONS[binarization_id]["transform"]
    except KeyError:
        log.error('Invalid binarization method')
        raise KeyError('Invalid binarization method')

    bin_matrices, additional_info = binarizer(matrix, axis=axis)
    
    if not isinstance(bin_matrices, List):
        bin_matrices = [bin_matrices]

    return bin_matrices, additional_info

def gt_val_to_gt_char(v):

    # assert (v >= 48 and v <= 57) or (v == -1)
    assert v >= -1

    if v == -1:
        # return chr(46) # '.'
        return 46
    else:
        return v+48

def recon_gt_mat_using_phasing_mat(
    allele_mat,
    phasing_mat,
    p
):
    pass

def reconstruct_genotype_matrix_using_phase_value(
    allele_mat, 
    phase_val, 
    p
):

    n_cols = allele_mat.shape[1]
    n_variants = allele_mat.shape[0]
    n_samples = n_cols//p
    phase_char = 47 # '/'

    max_val = allele_mat.max()
    char_per_genotype = len(str(max_val)) # In case #ALT > 9

    # ? number of genotypes{., 0, 1, ...} + number of phasing{\, |}  + number of separator {\t, \n}
    gt_mat_len = (n_cols*char_per_genotype + (p-1)*n_samples + n_samples) * n_variants

    n = 0
    curr_gt_val = 0
    # curr_sep
    variant_max_val = 0

    if phase_val == 0:
        phase_char = 124 # '|'

    gt_mat = [bytearray(b'\0' * gt_mat_len)]

    for i_variant in range(n_variants):

        variant_max_val = allele_mat[i_variant, :].max()

        # Use cheaper method using char
        if variant_max_val < 10:
            for i_sample in range(n_samples):

                curr_gt_val = allele_mat[i_variant, i_sample*p]
                
                # Handle case where GT for current sample is available
                if curr_gt_val == -2:
                    continue

                gt_mat[0][n] = gt_val_to_gt_char(curr_gt_val)
                n += 1

                for k in range(1, p):
                    curr_gt_val = allele_mat[i_variant, i_sample*p+k]
                    if curr_gt_val == -2:
                        break

                    gt_mat[0][n] = phase_char
                    gt_mat[0][n+1] = gt_val_to_gt_char(curr_gt_val)
                    n += 2
                
                # Add separator at the end of processing one sample
                if i_sample < n_samples-1:
                    gt_mat[0][n] = 9 # '\t'
                else:
                    gt_mat[0][n] = 10 # '\n'
                n += 1

        # Use more expensive method using string
        else:
            raise NotImplementedError()

    # length[0] = n
    length = n

    return gt_mat[0], length