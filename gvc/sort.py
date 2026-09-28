import time
import logging as log

import numpy as np
from scipy.spatial.distance import pdist, squareform
from numba import jit

import gvc.data_structures
from .utils import catchtime
from .dist import DIST_FUNC, comp_cost_mat
from .solver import SOLVERS
from .common import PERMUTATION_DTYPE

def _sort_matrix(
    bin_mat,
    dist_f_name=None,
    solver_name:str=None,
    solver_profile=0,
    sort_row=False,
    sort_col=False,
    transpose=False,
    **kwargs
):

    """
    Sort either row or column or both row and column of a matrix.
    Implements 4.3

    Parameters
    ----------
    bin_matrix : ndarray
        The binary matrix
    dist : string, optional, default: 'hrl'
        Distance function for cost matrix
    sort_row : bool, optional, default: False
        Sort rows of the matrix
    sort_col : bool, optional, default: False
        Sort columns of the matrix
    # transpose : bool, optional, default: False
    #     Transpose resulting sorted matrix
    time_limit : unsigned int, optional, default : None
        Limit the runtime of finding best tour.
        Set to None to disable time limit

    Returns
    -------
    sorted_matrix : ndarray (2d)
        Sorted matrix
    row_index : ndarray (1d)
        Indices stored in vector, used to reconstruct original matrix
        Value is None, if rows are not sorted
    col_index : ndarray (1d)
        Indices stored in vector, used to reconstruct original matrix
        Value is None, if columns are not sorted
    """
    
    assert solver_name in SOLVERS, "Invalid solver!"
    solver = SOLVERS[solver_name](solver_profile)
    
    if solver.req_dist_mat:
        assert dist_f_name in DIST_FUNC, "Invalid distance function!"

    if transpose:
        bin_mat = bin_mat.T
    else:
        bin_mat = bin_mat

    sorted_bin_mat = bin_mat.copy()

    log.debug("Matrix shape: {}".format(bin_mat.shape))
    
    if sort_row:
        log.info('Sort rows')
        
        if solver.req_dist_mat:
            with catchtime() as t:
                cost_mat = comp_cost_mat(bin_mat, dist_f_name)
            log.debug("Cost matrix time : {:.2f}s".format(t.time))
            
            with catchtime() as t:
                row_ids = solver.sort_rows(cost_mat)
            log.debug("Sort time: {:.2f}s".format(t.time))
            
        else:
            with catchtime() as t:
                row_ids = solver.sort_rows(bin_mat)
            log.debug("Total time: {:.2f}s".format(t.time))
        
        sorted_bin_mat = sorted_bin_mat[row_ids, :]
        row_ids = np.argsort(row_ids).astype(PERMUTATION_DTYPE)
        
    else:
        row_ids = None
    
    if sort_col:
        log.info('Sort columns')
  
        if solver.req_dist_mat:
            with catchtime() as t:
                cost_mat = comp_cost_mat(bin_mat.T, dist_f_name)
            log.debug("Cost matrix time : {:.2f}s".format(t.time))
            
            with catchtime() as t:
                col_ids = solver.sort_cols(cost_mat)
            log.debug("Sort time: {:.2f}s".format(t.time))
            
        else:
            with catchtime() as t:
                col_ids = solver.sort_cols(bin_mat)
            log.debug("Total time: {:.2f}s".format(t.time))

        sorted_bin_mat = sorted_bin_mat[:, col_ids]
        col_ids = np.argsort(col_ids).astype(PERMUTATION_DTYPE)
    else:
        col_ids = None

    return sorted_bin_mat, row_ids, col_ids

def sort(
    param_set:gvc.data_structures.ParameterSet,
    bin_allele_matrices:np.ndarray,
    phasing_matrix:np.ndarray,
    dist_f_name=None,
    solver_name:str=None,
    solver_profile=0,
):

    """
    Sort either row or column or both row and column of both matrices.
    Implements 4.3

    Parameters
    ----------
    param_set: gvc.data_structures.ParameterSet
        Parameter set
    bin_allele_matrices: np.ndarray, List of np.ndarray
        Binarized allele matrix
        List of ndarray when binarization using bit plane and concat_axis equals 2
    phasing_matrix: np.ndarray
        Phasing matrix
    dist: string, optional, default: 'hrl'
        Distance function for cost matrix
    time_limit : unsigned int, optional, default : None
        Limit the runtime of finding best tour.
        Set to None to disable time limit

    Returns
    -------
    data: dict

    """

    sorted_allele_matrices = []
    row_idx_allele_matrices = []
    col_idx_allele_matrices = []

    log.info('Sort variant matrices')
    for i in range(param_set.num_variants_flags):

        log.debug('Sort variant matrix #{}'.format(i))

        sorted_matrix, row_index, col_index = _sort_matrix(
            bin_allele_matrices[i],
            dist_f_name=dist_f_name,
            solver_name=solver_name,
            solver_profile=solver_profile,
            sort_row=param_set.sort_variants_row_flags[i],
            sort_col=param_set.sort_variants_col_flags[i],
            transpose=param_set.transpose_variants_mat_flags[i],
        )

        sorted_allele_matrices.append(sorted_matrix)
        row_idx_allele_matrices.append(row_index)
        col_idx_allele_matrices.append(col_index)

    if param_set.encode_phase_data:
        log.info('Sort phase matrix')
        sorted_phase_matrix, row_idx_phase_matrix, col_idx_phase_matrix = _sort_matrix(
            phasing_matrix,
            dist_f_name=dist_f_name,
            solver_name=solver_name,
            solver_profile=solver_profile,
            sort_row=param_set.sort_phases_row_flag,
            sort_col=param_set.sort_phases_col_flag,
            transpose=param_set.transpose_phase_mat_flag,
        )

    else:
        log.debug('Phase matrix is not sorted')
        sorted_phase_matrix = None
        row_idx_phase_matrix = None
        col_idx_phase_matrix = None

    return [
        sorted_allele_matrices,
        row_idx_allele_matrices,
        col_idx_allele_matrices,
        sorted_phase_matrix,
        row_idx_phase_matrix,
        col_idx_phase_matrix,
    ]