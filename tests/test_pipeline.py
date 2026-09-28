from io import BytesIO

import numpy as np
import pytest

from gvc.codec import MAT_CODECS
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.decoder import decode_encoded_variants
from gvc.encoder import run_core


def _array_encode(matrix):
    buffer = BytesIO()
    np.save(buffer, np.asarray(matrix), allow_pickle=False)
    return buffer.getvalue()


def _array_decode(payload):
    return np.load(BytesIO(payload), allow_pickle=False)


@pytest.fixture
def lossless_test_codec(monkeypatch):
    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", _array_encode)
    monkeypatch.setitem(codec, "decoder", _array_decode)


@pytest.mark.parametrize(
    "binarization_id,axis",
    [
        (BinarizationID.BIT_PLANE, 0),
        (BinarizationID.BIT_PLANE, 2),
        (BinarizationID.ROW_BIN_SPLIT, 0),
    ],
)
@pytest.mark.parametrize(
    "sort_rows,sort_cols",
    [
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ],
)
def test_core_pipeline_roundtrip(
    lossless_test_codec,
    binarization_id,
    axis,
    sort_rows,
    sort_cols,
):
    allele_matrix = np.array(
        [
            [0, 1, 1, 1, 2, 0],
            [3, 1, 0, 2, 1, 1],
            [1, 0, 2, 2, 0, 3],
            [0, 0, 1, 2, 3, 1],
        ],
        dtype=np.uint8,
    )
    phase_matrix = np.array(
        [
            [False, True, False],
            [True, False, True],
            [False, False, True],
            [True, True, False],
        ],
        dtype=bool,
    )

    raw_block = (allele_matrix.copy(), phase_matrix.copy(), 2, None, None)
    ps_params = [
        binarization_id,
        CodecID.JBIG1,
        axis,
        sort_rows,
        sort_cols,
        False,
    ]
    tsp_params = ["ham", "nn", 0]

    block, parameter_set = run_core(raw_block, ps_params, tsp_params)
    restored_alleles, restored_phase = decode_encoded_variants(
        parameter_set,
        block.block_payload,
        ret_gt=False,
    )

    np.testing.assert_array_equal(restored_alleles, allele_matrix)
    np.testing.assert_array_equal(restored_phase, phase_matrix)
