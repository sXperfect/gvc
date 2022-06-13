
import numpy as np
import logging as log

from .data_structures import AccessUnit, ParameterSet
from .binarization import binarization_str_to_flag
from .codec import CODEC_STR2ID

signed_allele_dtype = np.int32
allele_dtype = np.uint32

phasing_dtype = np.bool

max_val_dtype = allele_dtype
bin_dtype = np.bool

def create_parameter_set(
    any_missing,
    not_available,
    p:int,
    phasing_matrix,
    additional_info,
    binarization,
    codec_name,
    axis,
    sort_rows,
    sort_cols,
    transpose=False,
    parameter_set_id=0
):

    # Convert from string to int
    binarization_flag = binarization_str_to_flag[binarization]

    codec_id = CODEC_STR2ID[codec_name]

    # Binarization using bit plane
    if binarization_flag == 0:
        num_bin_mat = int(additional_info)

        if axis not in (0,1,2):
            log.error('Invalid value for concat_axis:{}'.format(axis))
            raise ValueError('Invalid value for concat_axis:{}'.format(axis))

        if axis in (0, 1):
            num_variants_flags = 1

        elif axis == 2:
            num_variants_flags = num_bin_mat

    # Binarization by row splitting
    elif binarization_flag in (1, 2, 3):
        num_bin_mat = 1
        num_variants_flags = num_bin_mat

    else:
        log.error('Invalid binarization_flag: {}'.format(binarization_flag))
        raise ValueError('Invalid binarization_flag: {}'.format(binarization_flag))

    # Handle case where phasing matrix can be represented by a single value
    if p == 1 or np.all(~phasing_matrix) or np.all(phasing_matrix):
        if phasing_matrix is None:
            phasing_value = 0
        else:
            # First element is sufficient to represent the matrix
            phasing_value = phasing_matrix[0,0]

        param_set = ParameterSet(
            parameter_set_id,
            any_missing,
            not_available,
            p,
            binarization_flag,
            num_bin_mat,
            axis,      # concat_axis
            [sort_rows for _ in range(num_variants_flags)],  # sort_variants_row_flags
            [sort_cols for _ in range(num_variants_flags)],  # sort_variants_col_flags
            [transpose for _ in range(num_variants_flags)],  # transpose_variants_mat_flags
            [codec_id for _ in range(num_variants_flags)],
            False,  # Do not encode phase data
            phase_value=phasing_value
        )
    else:
        param_set = ParameterSet(
            parameter_set_id,
            any_missing,
            not_available,
            p,
            binarization_flag,
            num_bin_mat,
            axis,      # concat_axis
            [sort_rows for _ in range(num_variants_flags)],  # sort_variants_row_flags
            [sort_cols for _ in range(num_variants_flags)],  # sort_variants_col_flags
            [transpose for _ in range(num_variants_flags)],  # transpose_variants_mat_flags
            [codec_id for _ in range(num_variants_flags)],
            True,  # Do not encode phase data
            sort_phases_row_flag=sort_rows,
            sort_phases_col_flag=sort_cols,
            transpose_phase_mat_flag=transpose,
            phase_coder_ids=codec_id
        )

    return param_set

def store_access_unit(
    f,
    access_unit_id,
    parameter_set,
    blocks,
):
    acc_unit = AccessUnit.from_blocks(
        access_unit_id, 
        parameter_set.parameter_set_id,
        blocks
    )

    f.write(acc_unit.to_bytes())