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
    tensor = np.asarray(tensor)
    if tensor.ndim != 3:
        raise ValueError("tensor must be three-dimensional")
    if tensor.shape[1] <= 0:
        raise ValueError("tensor sample dimension must be positive")

    list_matrix = np.split(tensor, tensor.shape[1], axis=1)
    matrix = np.concatenate(list_matrix, axis=2).squeeze(axis=1)
    if matrix.dtype != tensor.dtype:
        raise RuntimeError("tensor-to-matrix conversion changed dtype")
    return matrix

def simd_tensor_to_txt(allele_tensor, phasing_tensor):
    """Render genotype tensors without single-character allele assumptions."""
    allele_tensor = np.asarray(allele_tensor)
    phasing_tensor = np.asarray(phasing_tensor)

    if allele_tensor.ndim != 3:
        raise ValueError("allele tensor must be three-dimensional")
    if phasing_tensor.ndim != 3:
        raise ValueError("phasing tensor must be three-dimensional")

    n_variants, n_samples, ploidy = allele_tensor.shape
    expected_phase_shape = (n_variants, n_samples, max(0, ploidy - 1))
    if phasing_tensor.shape != expected_phase_shape:
        raise ValueError(
            "phasing tensor shape does not match allele tensor/ploidy"
        )

    rows = []
    for i_variant in range(n_variants):
        genotypes = []
        for i_sample in range(n_samples):
            alleles = allele_tensor[i_variant, i_sample]
            if int(alleles[0]) == -2:
                genotypes.append("")
                continue

            text = allele_val2str(int(alleles[0]))
            for k in range(1, ploidy):
                allele = int(alleles[k])
                if allele == -2:
                    break
                phase = int(phasing_tensor[i_variant, i_sample, k - 1])
                if phase not in (0, 1):
                    raise ValueError("phasing value must be 0 or 1")
                text += PHASING_VAL2CHAR[phase]
                text += allele_val2str(allele)
            genotypes.append(text)
        rows.append("\t".join(genotypes))
    return "\n".join(rows) + ("\n" if rows else "")

def matrix_to_tensor(matrix, num_matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if not isinstance(num_matrix, int) or isinstance(num_matrix, bool) or num_matrix <= 0:
        raise ValueError("num_matrix must be a positive integer")
    if matrix.shape[1] % num_matrix:
        raise ValueError("matrix column count must be divisible by num_matrix")

    list_matrix = np.split(
        np.expand_dims(matrix, axis=1),
        matrix.shape[1] // num_matrix,
        axis=2,
    )

    return np.concatenate(list_matrix, axis=1)

PHASING_VAL2CHAR = ['|', '/']
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
    