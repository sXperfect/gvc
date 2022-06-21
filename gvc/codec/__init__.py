import logging as log
import typing as t
import numpy as np

from ..data_structures.consts import CodecID
from ..data_structures import RowColIds, ParameterSet, VectorAMax

from . import jbigkit #? Import jbigkit

#? If a new codec is added, please update data_structure.consts too
MAT_CODECS = {
    CodecID.JBIG1 : {
        "name": "jbig",
        "encoder": jbigkit.encode, #? Add encode function
        "decoder": jbigkit.decode, #? Add decode function
    }
}

AVAIL_CODECS = [v['name'] for v in MAT_CODECS.values()]
CODEC_STR2ID = {v["name"]:k for k, v in MAT_CODECS.items()}

def decode_permutation(payload, num_entries):
    permutation = RowColIds.from_bytes(payload, num_entries).ids
    return permutation

def encode_permutation(permutation):
    payload = RowColIds(permutation).to_bitio().to_bytes(align=True)
    return payload

def decode_vector(payload):
    vect = VectorAMax.from_bytes(payload).vector
    return vect

def encode_vector(vect):
    payload = VectorAMax(vect).to_bitio().to_bytes(align=True)
    return payload

def decode(
    bin_mat_payload:bytes, 
    row_ids_payload:bytes, 
    col_ids_payload:bytes, 
    decoder_id:int,
    unsort=True
):
    try:
        bin_mat_payload = bin_mat_payload.read()
    except:
        pass

    decode_f = MAT_CODECS[decoder_id]["decoder"]
    bin_mat = decode_f(bin_mat_payload, decoder_id)
    nrows, ncols = bin_mat.shape
    
    if row_ids_payload is not None:
        try:
            row_ids_payload = row_ids_payload.read()
        except:
            pass
        
        row_ids = decode_permutation(row_ids_payload, nrows)
    else:
        row_ids = None
    
    if col_ids_payload:
        try:
            col_ids_payload = col_ids_payload.read()
        except:
            pass

        col_ids = decode_permutation(col_ids_payload, ncols)
    else:
        col_ids = None
        
    if unsort:
        return bin_mat
    else:
        return bin_mat, row_ids, col_ids

def encode(
    param_set:ParameterSet,
    additional_info,
    sorted_allele_matrices:t.List[np.ndarray],
    row_idx_allele_matrices:t.List[np.ndarray],
    col_idx_allele_matrices:t.List[np.ndarray],
    sorted_phase_mat:np.ndarray,
    row_idx_phase_matrix:np.ndarray,
    col_idx_phase_matrix:np.ndarray,
):


    """

    Parameters
    ----------
    param_set:gvc.data_structures.ParameterSet,
        The parameter set for current payload.
    additional_info: #TODO: Add desc
        #TODO: Add desc
    sorted_allele_matrices: List[np.ndarray],
        List of sorted allele matrices.
    row_idx_allele_matrices: List[np.ndarray],
        List of row index, entriy of the list is None if no sorting was done
    col_idx_allele_matrices: List[np.ndarray],
        List of col index, entriy of the list is None if no sorting was done
    sorted_phase_matrix: np.ndarray,
        The sorted phase matrix
    row_idx_phase_matrix: np.ndarray,
        Row index of phase matrix, if phase matrix is to be encoded and the rows are sorted
    col_idx_phase_matrix: np.ndarray,
        Col index of phase matrix, if phase matrix is to be encoded and the columns are sorted.

    Returns
    -------

    """
    
    sorted_allele_mat_payloads = []
    row_idx_allele_mat_payloads = []
    col_idx_allele_mat_payloads = []
    for i, matrix in enumerate(sorted_allele_matrices):

        log.info('Encode matrix number {}'.format(i))
        encoder_id = param_set.variants_coder_ids[i]
        encoder_f = MAT_CODECS[encoder_id]["encoder"]
        matrix_bytes = encoder_f(matrix)

        sorted_allele_mat_payloads.append(matrix_bytes)

        row_ids_payload = None
        if param_set.sort_variants_row_flags[i]:
            log.info('Encode row ids number {}'.format(i))
            row_ids_payload = encode_permutation(row_idx_allele_matrices[i])
            

        row_idx_allele_mat_payloads.append(row_ids_payload)

        col_ids_payload = None
        if param_set.sort_variants_col_flags[i]:
            log.info('Encode column ids number {}'.format(i))
            col_ids_payload  = encode_permutation(col_idx_allele_matrices[i])

        col_idx_allele_mat_payloads.append(col_ids_payload)
        
    if param_set.binarization_flag in (0, 3):
        variants_amax_payload = None
    elif param_set.binarization_flag in (1,2):
        variants_amax_payload = encode_vector(additional_info)
    else:
        raise ValueError('Invalid binarization flag')

    # Initialize default values
    sorted_phase_mat_payload = None
    row_ids_phase_mat_payload = None
    col_idx_phase_mat_payload = None

    if param_set.encode_phase_data:
        encoder_id = param_set.phase_coder_ids
        encoder_f = MAT_CODECS[encoder_id]["encoder"]
        sorted_phase_mat_payload = encoder_f(sorted_phase_mat)

        if param_set.sort_phases_row_flag:
            log.info('Encode row ids of phase matrix')
            row_ids_phase_mat_payload = encode_permutation(row_idx_phase_matrix)

        if param_set.sort_phases_col_flag:
            log.info('Encode column ids of phase matrix')
            col_idx_phase_mat_payload = encode_permutation(col_idx_phase_matrix)

    return [
        sorted_allele_mat_payloads,
        row_idx_allele_mat_payloads,
        col_idx_allele_mat_payloads,
        variants_amax_payload,
        sorted_phase_mat_payload,
        row_ids_phase_mat_payload,
        col_idx_phase_mat_payload,
    ]
    
def decode(
    bin_mat_payload:bytes, 
    row_ids_payload:bytes, 
    col_ids_payload:bytes, 
    coder_id:int,
    unsort=True
):
    """

    Parameters
    ----------
    bin_mat_payload : RandomAccessHandler or bytes
    row_ids_payload : RandomAccessHandler or bytes, 
    col_ids_payload : RandomAccessHandler or bytes, 
    coder_id:int
    Returns
    -------

    """
    try:
        bin_mat_payload = bin_mat_payload.read()
    except AttributeError:
        pass

    decoder_f = MAT_CODECS[coder_id]['decoder']
    bin_mat = decoder_f(bin_mat_payload)
    nrows, ncols = bin_mat.shape

    if row_ids_payload is not None:
        try:
            row_ids_payload = row_ids_payload.read()
        except AttributeError:
            pass
        row_ids = decode_permutation(row_ids_payload, nrows)

        log.info('Unsort rows')

        if unsort:
            bin_mat = bin_mat[row_ids, :]
    else:
        row_ids = None

    if col_ids_payload is not None:
        try:
            col_ids_payload = col_ids_payload.read()
        except AttributeError:
            pass

        col_ids = decode_permutation(col_ids_payload, ncols)

        log.info('Unsort columns')

        if unsort:
            bin_mat = bin_mat[:, col_ids]
    else:
        col_ids = None

    if unsort:
        return bin_mat
    else:
        return bin_mat, row_ids, col_ids