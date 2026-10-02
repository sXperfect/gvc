"""Matrix ordering for the GVC transform pipeline."""

import logging as log

import numpy as np

import gvc.data_structures
from .common import PERMUTATION_DTYPE
from .dist import DIST_FUNC, comp_cost_mat
from .solver import SOLVERS
from .utils import catchtime


def _sort_matrix(
    bin_mat,
    dist_f_name="ham",
    solver_name="nn",
    solver_profile=0,
    sort_row=False,
    sort_col=False,
    transpose=False,
    **kwargs
):
    if solver_name not in SOLVERS:
        raise ValueError("Invalid solver: {}".format(solver_name))

    solver = SOLVERS[solver_name](solver_profile)
    if solver.req_dist_mat and dist_f_name not in DIST_FUNC:
        raise ValueError("Invalid distance function: {}".format(dist_f_name))

    working = bin_mat.T if transpose else bin_mat
    sorted_bin_mat = working.copy()
    log.debug("Matrix shape: %s", working.shape)

    max_perm_entries = np.iinfo(PERMUTATION_DTYPE).max
    if sort_row and working.shape[0] > max_perm_entries:
        raise ValueError(
            "row sorting exceeds permutation format limit of {}".format(
                max_perm_entries
            )
        )
    if sort_col and working.shape[1] > max_perm_entries:
        raise ValueError(
            "column sorting exceeds permutation format limit of {}".format(
                max_perm_entries
            )
        )

    if sort_row:
        if solver.req_dist_mat:
            with catchtime() as timer:
                cost_mat = comp_cost_mat(working, dist_f_name)
            log.debug("Row cost matrix time: %.2fs", timer.time)
            row_ids = solver.sort_rows(cost_mat)
        else:
            row_ids = solver.sort_rows(working)
        sorted_bin_mat = sorted_bin_mat[row_ids, :]
        row_ids = np.argsort(row_ids).astype(PERMUTATION_DTYPE)
    else:
        row_ids = None

    if sort_col:
        if solver.req_dist_mat:
            with catchtime() as timer:
                cost_mat = comp_cost_mat(working.T, dist_f_name)
            log.debug("Column cost matrix time: %.2fs", timer.time)
            col_ids = solver.sort_cols(cost_mat)
        else:
            col_ids = solver.sort_cols(working)
        sorted_bin_mat = sorted_bin_mat[:, col_ids]
        col_ids = np.argsort(col_ids).astype(PERMUTATION_DTYPE)
    else:
        col_ids = None

    return sorted_bin_mat, row_ids, col_ids


def sort(
    param_set: gvc.data_structures.ParameterSet,
    bin_allele_matrices,
    phasing_matrix,
    dist_f_name="ham",
    solver_name="nn",
    solver_profile=0,
):
    sorted_allele_matrices = []
    row_idx_allele_matrices = []
    col_idx_allele_matrices = []

    for i in range(param_set.num_variants_flags):
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
