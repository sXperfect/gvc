from typing import List
import logging as log
import time

import numpy as np
import itertools as it

import gvc.common

from .data_structures import consts
from . import cdebinarize

_phasing_dict = {
    '|' : 0,
    '/' : 1
}

_phasing_val_to_str = ['|', '/']

def _gt_code_to_int(gt_code):
    try:
        return gvc.common.allele_dtype(gt_code)
    except ValueError:
        if gt_code == ".":
            return gvc.common.signed_allele_dtype(-1)
        else:
            raise ValueError()

def _int_to_gt_code(val):

    if val >= 0:
        return str(val)
    else:
        return "."

def _matrix_to_tensor(matrix, num_matrix):

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

def simd_tensor_to_txt(allele_tensor, phasing_tensor):
    m, n, p = allele_tensor.shape

    max_val = allele_tensor.max()
    avail_allele_vals = '.' + "".join(np.arange(max_val+1).astype(str))
    avail_phase_val = '|/'

    if p == 1:
        codebook = avail_allele_vals
    elif p > 1:
        codebook = ["".join(l) for l in it.product(avail_allele_vals, avail_phase_val, repeat=p-1)]
        codebook = ["".join(l) for l in it.product(codebook, avail_allele_vals)]

    def_val_idx = 0
    for k in range(1, p):
        def_val_idx += len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k)
        # def_val_idx += 0*len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k-1)
    def_val_idx += 1

    allele_tensor += 1

    flag_mask = allele_tensor == -1
    complete_allele_mask = (np.logical_not(flag_mask)).all(axis=2)

    # gt_matrix_str = np.full((m,n), codebook[def_val_idx])
    gt_matrix_str = np.empty((m,n), dtype='<U{}'.format(2*p))

    index_mat = np.zeros(gt_matrix_str.shape, dtype=np.int)
    for k in range(1, p):
        allele_factor = len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k)
        index_mat[complete_allele_mask] += allele_tensor[complete_allele_mask, k-1] * allele_factor
        phase_factor = len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k-1)
        index_mat[complete_allele_mask] += phasing_tensor[complete_allele_mask, k-1] * phase_factor

    index_mat[complete_allele_mask] += allele_tensor[complete_allele_mask, -1]
    gt_matrix_str[complete_allele_mask] = np.array(codebook)[index_mat[complete_allele_mask]]

    if (flag_mask).any():
        for i, j, k in np.argwhere(flag_mask):
            sample_val = avail_allele_vals[allele_tensor[i, j, 0]]
            for l in range(1, k):
                sample_val = avail_allele_vals[allele_tensor[i, j, l]] + avail_phase_val[phasing_tensor[i, j, l-1]]

            gt_matrix_str[i, j] = sample_val

    delim_mat_str = np.full((m,n,1), '\t')
    delim_mat_str[:, -1] = '\n'
    gt_matrix_str = np.concatenate((np.expand_dims(gt_matrix_str,2), delim_mat_str), axis=2)
    gt_matrix_str = _tensor_to_matrix(gt_matrix_str)
    
    return "".join(gt_matrix_str.reshape(-1))

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

    allele_tensor = _matrix_to_tensor(allele_matrix, p)

    if p >= 1:
        if isinstance(phasing_matrix, int):
            allele_tensor_shape = allele_tensor.shape

            phasing_tensor = np.zeros((*allele_tensor_shape[:-1], p-1), dtype=np.bool)
            phasing_tensor[:, :, :] = phasing_matrix

        else:
            phasing_tensor = _matrix_to_tensor(phasing_matrix, p-1)

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
        any_missing = True

        new_max_val = np.max(allele_matrix) + 1
        allele_matrix[dot_mask] = new_max_val
    else:
        any_missing = False

    if np.any(na_mask):
        not_available = True

        new_max_val = np.max(allele_matrix) + 1
        allele_matrix[na_mask] = new_max_val
    else:
        not_available = False

    return allele_matrix.astype(np.uint8), any_missing, not_available

def undo_adaptive_max_value(allele_array, any_missing_flag, not_available_flag):
    allele_array = allele_array.astype(gvc.common.signed_allele_dtype)

    if not_available_flag:
        true_max_val = np.max(allele_array)

        allele_array[allele_array == true_max_val] = -2

    if any_missing_flag:
        true_max_val = np.max(allele_array)

        allele_array[allele_array == true_max_val] = -1
    
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

    allele_tensor = np.zeros([m, n, p], dtype=gvc.common.signed_allele_dtype)
    assert(np.issubdtype(allele_tensor.dtype, np.signedinteger))

    if p-1 > 0:
        phasing_tensor = np.zeros([m, n, (p - 1)], dtype=gvc.common.phasing_dtype)
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
                            -2*np.ones([m, n, p_new - p], dtype=gvc.common.signed_allele_dtype)
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


def _binarize_using_bit_plane(matrix, axis=None, **kwargs):
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

    bit_depth = np.ceil(np.log2(matrix.max() + 1)).astype(gvc.common.allele_dtype)

    for i_bit in range(bit_depth):
        bit_tensor = np.bitwise_and(matrix, int(2**i_bit)).astype(gvc.common.bin_dtype)
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

def _debinarize_using_bit_plane(bin_matrices, bit_depth, axis):

    if axis < 2:
        assert len(bin_matrices) == 1

        bit_planes = np.split(bin_matrices[0], bit_depth, axis)

    elif axis == 2 :
        bit_planes = bin_matrices

    else:
        log.error('Invalid axis: {}'.format(axis))
        raise gvc.errors.GvcError()

    assert bit_depth == len(bit_planes)

    matrix = bit_planes[0].astype(gvc.common.allele_dtype)
    for i in range(1, bit_depth):
        matrix += 2**i * bit_planes[i].astype(gvc.common.allele_dtype)

    return matrix

def _binarize_by_row_splitting(matrix, **kwargs):
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
    max_val_per_row : ndarray (1d)
        Maximum value on each row 
    """

    log.debug('binarizing allele matrix by row splitting')

    nrow = matrix.shape[0]

    max_val_per_row = np.max(matrix, axis=1).astype(gvc.common.max_val_dtype)
    # Force max value 0 to 1
    max_val_per_row[max_val_per_row == 0] = 1

    # consider the case when a line has all 0s
    log.info("Greatest value in the block: {}".format(max_val_per_row.max()))
    bin_mat_shape = (np.sum(max_val_per_row).astype(int), matrix.shape[1])
    bin_mat = np.zeros(bin_mat_shape, dtype=gvc.common.bin_dtype)

    r = 0  # row id in bin_mat
    inBinSubblock = False
    binSubblockStart = 0
    for i in range(nrow):  # iterate over rows in matrix

        if max_val_per_row[i] <= 1:  # if row of matrix is binary: copy
            if not inBinSubblock:
                inBinSubblock = True
                binSubblockStart = i

        else:

            # end of binary subblock
            if inBinSubblock:
                inBinSubblock = False
                binSubblockSize = i - binSubblockStart
                bin_mat[r:r + binSubblockSize, :] = matrix[i - binSubblockSize:i, :]
                r += binSubblockSize

            # if row of matrix is non-binary, binarize only row and append max_val_per_row[i] rows to bin_mat
            for rPrime in range(max_val_per_row[i]):
                bin_mat[r] = (matrix[i] == (rPrime+1))
                r += 1

    if inBinSubblock:
        binSubblockSize = nrow - binSubblockStart
        bin_mat[r:r + binSubblockSize, :] = matrix[nrow - binSubblockSize:nrow, :]

    return bin_mat, max_val_per_row

def _debinarize_by_row_splitting(bin_matrices, amax, *args):

    try:
        if isinstance(bin_matrices, List):
            bin_mat = bin_matrices[0]
        elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 3:
            bin_mat = bin_matrices[0, :, :]
        elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 2:
            bin_mat = bin_matrices
        else:
            raise TypeError
    except TypeError:
        log.info("Invalid datat type:{}".format(type(bin_matrices)))
        raise ValueError

    nrows = amax.shape[0]
    __, ncols = bin_mat.shape

    matrix = np.zeros((nrows, ncols), dtype=gvc.common.allele_dtype)

    irow_bin_mat = 0
    for i in range(nrows):
        for j in range(amax[i]):
            # matrix[i, :] <<= 1
            matrix[i, :] += ((j+1) * bin_mat[irow_bin_mat+j, :]).astype(gvc.common.allele_dtype)

        irow_bin_mat += amax[i]

    return matrix

def _binarize_by_rc_bin_split(matrix, axis=0, **kwargs):
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

    if axis == 1:
        matrix = matrix.T

    nrow = matrix.shape[0]

    bitlen_vect = np.max(matrix, axis=1).astype(np.uint)
    # Force max value 0 to 1
    bitlen_vect[bitlen_vect == 0] = 1

    # consider the case when a line has all 0s
    log.info("Greatest value in the block: {}".format(bitlen_vect.max()))

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

    if axis == 1:
        bin_mat = bin_mat.T

    return bin_mat, bitlen_vect


def _debinarize_by_rc_bin_split(bin_matrices, bitlen_vect, **kwargs):

    try:
        if isinstance(bin_matrices, List):
            bin_mat = bin_matrices[0]
        elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 3:
            bin_mat = bin_matrices[0, :, :]
        elif isinstance(bin_matrices, np.ndarray) and bin_matrices.ndim == 2:
            bin_mat = bin_matrices
        else:
            raise TypeError
    except TypeError:
        log.info("Invalid datat type:{}".format(type(bin_matrices)))
        raise ValueError

    # if axis == 1:
    #     bin_mat = bin_mat.T

    matrix = cdebinarize.debin_rc_bin_split(bin_mat, bitlen_vect.astype(np.uint8))
    # matrix = gvc.debinarize.debin_rc_bin_split(bin_mat, bitlen_vect.astype(np.uint8))

    # if axis == 1:
    #     matrix = matrix.T

    return matrix

def _binarize_rc_bin_split_v2(matrix, axis=0, **kwargs):
    """
    Binarize a matrix by row column binary splitting.

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

    if axis == 1:
        matrix = matrix.T

    nrow = matrix.shape[0]

    bitlen_vect = np.max(matrix, axis=1).astype(np.uint)
    # Force max value 0 to 1
    bitlen_vect[bitlen_vect == 0] = 1

    # consider the case when a line has all 0s
    log.info("Greatest value in the block: {}".format(bitlen_vect.max()))

    bitlen_vect = np.ceil(np.log2(bitlen_vect + 1)).astype(np.uint16)
    bin_mat_nrows = np.sum(bitlen_vect).astype(np.uint)
    bin_mat_shape = (bin_mat_nrows, matrix.shape[1]+1)
    bin_mat = np.zeros(bin_mat_shape, dtype=np.bool)

    i_b = 0  # row id in bin_mat
    for i in range(nrow):  # iterate over rows in matrix

        bitlen = bitlen_vect[i]

        # bin_mat[i_b, 0] = True
        for i_bit in range(bitlen):
            bin_mat[i_b + i_bit, 1:] = np.bitwise_and(matrix[i, :], 1<<i_bit).astype(np.bool)
        
        bin_mat[i_b+i_bit, 0] = True

        i_b += bitlen

    if axis == 1:
        bin_mat = bin_mat.T

    return bin_mat, None

avail_binarization = {
    'bit_plane' : _binarize_using_bit_plane,
    'row_split' : _binarize_by_row_splitting,
    'rc_bin_split' : _binarize_by_rc_bin_split,
    'rc_bin_split_2' : _binarize_rc_bin_split_v2,
}

BIT_PLANE_FLAG = 0
ROW_SPLIT_FLAG = 1
RC_BIN_SPLIT_FLAG = 2

avail_debinarization = {
    0 : _debinarize_using_bit_plane,
    1 : _debinarize_by_row_splitting,
    2 : _debinarize_by_rc_bin_split,
    3 : _debinarize_by_rc_bin_split,
}

binarization_str_to_flag = {
    'bit_plane' : 0,
    'row_split' : 1,
    'rc_bin_split' : 2,
    'rc_bin_split_2' : 3
}

def debinarize_bin_matrices(
    bin_matrices,
    additional_info,
    axis,
    binarization_flag:int
):
    try:
        debinarizer = avail_debinarization[binarization_flag]
    except IndexError:
        log.error('Invalid binarization method')
        raise KeyError('Invalid binarization method')

    return debinarizer(bin_matrices, additional_info, axis=axis)

def binarize_allele_matrix(
    matrix, 
    binarizer=None, 
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
        binarization_method = avail_binarization[binarizer]
    except KeyError:
        log.error('Invalid binarization method')
        raise KeyError('Invalid binarization method')

    bin_matrices, additional_info = binarization_method(matrix, axis=axis)
    
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

    # gt_mat[0] = <char *> malloc((gt_mat_len + 1) * sizeof(char))
    gt_mat = [bytearray(b'\0' * gt_mat_len)]
    # if not gt_mat[0]:
    #     raise MemoryError()

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