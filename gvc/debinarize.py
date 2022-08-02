import itertools as it
import numpy as np

def tensor_to_matrix(tensor):
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

def simd_tensor_to_txt(allele_tensor, phasing_tensor):
    m, n, p = allele_tensor.shape

    max_val = allele_tensor.max()
    avail_allele_vals = '.' + "".join(np.arange(max_val+1).astype(str))
    avail_phase_val = '/|'

    if p == 1:
        codebook = avail_allele_vals
    elif p > 1:
        codebook = ["".join(l) for l in it.product(avail_allele_vals, avail_phase_val, repeat=p-1)]
        codebook = ["".join(l) for l in it.product(codebook, avail_allele_vals)]

    def_val_idx = 0
    for k in range(1, p):
        def_val_idx += len(avail_allele_vals) ** (p-k) * len(avail_phase_val) ** (p-k)
    def_val_idx += 1

    allele_tensor += 1

    flag_mask = allele_tensor == -1
    complete_allele_mask = (np.logical_not(flag_mask)).all(axis=2)

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
    gt_matrix_str = tensor_to_matrix(gt_matrix_str)
    
    return "".join(gt_matrix_str.reshape(-1))

def matrix_to_tensor(matrix, num_matrix):

    list_matrix = np.split(
        np.expand_dims(matrix, axis=1), 
        matrix.shape[1]//num_matrix, 
        axis=2
    )

    return np.concatenate(list_matrix, axis=1)

PHASING_VAL2CHAR = ['/', '|']
ALLELE_VAL2CHAR = np.arange(18).astype(object)
ALLELE_VAL2CHAR[-2] = ''
ALLELE_VAL2CHAR[-1] = '.'

def debin_rc_bin_split(bin_mat, bitlen_vect):
    
    bitlen_vect = bitlen_vect.astype(np.uint8)
    
    nrows = len(bitlen_vect)
    ncols = bin_mat.shape[1]
    
    mat = np.zeros((nrows, ncols), dtype=np.uint8)
    
    irow_bin_mat = 0
    for i in range(nrows):
        for j in range(bitlen_vect[i]):
            
            int_row = bin_mat[irow_bin_mat+j, :].astype(np.uint8) << j
            mat[i, :] |= int_row

        irow_bin_mat += bitlen_vect[i]
        
    return mat

def allele_val2str(v):
    if v >= 0:
        return str(v)
    elif v == -1:
        return '.'
    else:
        return ''
    
def phasing_val2str(v):
    return PHASING_VAL2CHAR[v]
    
def my_func(v):
    
    return ''.join(v)
    
vectorized_allele_val2str = np.vectorize(allele_val2str)
vectorized_phasing_val2str = np.vectorize(phasing_val2str)

def recon_gt_mat_with_phase_val(allele_mat, phasing_val, p):
    
    # allele_min_val = allele_mat.min()
    # allele_max_val = allele_mat.max()
    
    allele_tensor = matrix_to_tensor(allele_mat, p)
    #TODO: Assume p is always greater than 1
    phasing_tensor = np.full(
        [*allele_tensor.shape[0:2], p-1],
        phasing_val
    )

    out = simd_tensor_to_txt(allele_tensor, phasing_tensor)
    
    return out
    
    # phasing_char = PHASING_VAL2CHAR[phasing_val]
    
    # # allele_str_mat = vectorized_allele_val2str(allele_mat)
    # allele_str_mat = ALLELE_VAL2CHAR[allele_mat]
    # allele_str_tensor = matrix_to_tensor(allele_str_mat, p)
    
    # for i_p in range(p-1):
    #     idx = i_p*p+1
        
    #     allele_str_tensor = np.insert(allele_str_tensor, idx, phasing_char, axis=2)
    
    # gt_mat = np.apply_along_axis(my_func, 2, allele_str_tensor)

    # return gt_mat

def recon_gt_mat_with_phase_mat(allele_mat, phasing_mat, p):
    
    allele_tensor = matrix_to_tensor(allele_mat, p)
    phasing_tensor = matrix_to_tensor(phasing_mat, p-1)
    
    out = simd_tensor_to_txt(allele_tensor, phasing_tensor)
    return out
    
    # allele_str_mat = vectorized_allele_val2str(allele_mat)
    # allele_str_tensor = matrix_to_tensor(allele_str_mat, p)
    
    # if p > 1:
    #     phasing_str_mat = vectorized_phasing_val2str(phasing_mat)
        
    #     if p > 2:
    #         phasing_str_tensor = matrix_to_tensor(phasing_str_mat, p-1)
    #     else:
    #         phasing_str_tensor = np.expand_dims(phasing_str_mat, -1)
            
    #     for i_p in range(p-1):
    #         idx = i_p*p+1
            
    #         allele_str_tensor = np.insert(allele_str_tensor, idx, phasing_str_tensor[:, :, i_p], axis=2)
            
    # gt_mat = np.apply_along_axis(my_func, 2, allele_str_tensor)
    
    # return gt_mat
    