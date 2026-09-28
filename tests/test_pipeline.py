from io import BytesIO

import numpy as np
import pytest

from gvc import binarization
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



def _make_raw_block(seed, ploidy, phase_mode):
    rng = np.random.default_rng(seed)
    n_variants = 5
    n_samples = 4

    source = rng.integers(
        0,
        5,
        size=(n_variants, n_samples * ploidy),
        dtype=np.int8,
    )
    # Exercise adaptive missing and not-available values without making a row
    # entirely unavailable.
    source[0, 0] = -1
    if source.shape[1] > 1:
        source[1, 1] = -2

    encoded, missing_rep, na_rep = binarization.adaptive_max_value(source)

    if ploidy == 1:
        phases = np.empty((n_variants, 0), dtype=bool)
    else:
        phase_cols = n_samples * (ploidy - 1)
        if phase_mode == "phased":
            phases = np.zeros((n_variants, phase_cols), dtype=bool)
        elif phase_mode == "unphased":
            phases = np.ones((n_variants, phase_cols), dtype=bool)
        else:
            phases = rng.integers(
                0,
                2,
                size=(n_variants, phase_cols),
                dtype=np.uint8,
            ).astype(bool)

    return source, (encoded, phases, ploidy, missing_rep, na_rep)


@pytest.mark.parametrize("ploidy", [1, 2, 3, 4])
@pytest.mark.parametrize("phase_mode", ["phased", "unphased", "mixed"])
@pytest.mark.parametrize(
    "binarization_id,axis,sort_rows,sort_cols,transpose",
    [
        (BinarizationID.BIT_PLANE, 0, False, False, False),
        (BinarizationID.BIT_PLANE, 1, True, False, False),
        (BinarizationID.BIT_PLANE, 2, False, True, False),
        (BinarizationID.BIT_PLANE, 2, True, True, True),
        (BinarizationID.ROW_BIN_SPLIT, 0, False, False, False),
        (BinarizationID.ROW_BIN_SPLIT, 0, True, True, True),
    ],
)
def test_generated_pipeline_roundtrip_across_v1_parameter_space(
    lossless_test_codec,
    ploidy,
    phase_mode,
    binarization_id,
    axis,
    sort_rows,
    sort_cols,
    transpose,
):
    source, raw_block = _make_raw_block(
        seed=3800 + ploidy * 10 + len(phase_mode),
        ploidy=ploidy,
        phase_mode=phase_mode,
    )

    block, parameter_set = run_core(
        raw_block,
        [
            binarization_id,
            CodecID.JBIG1,
            axis,
            sort_rows,
            sort_cols,
            transpose,
        ],
        ["ham", "nn", 0],
    )

    restored_alleles, restored_phase = decode_encoded_variants(
        parameter_set,
        block.block_payload,
        ret_gt=False,
    )
    np.testing.assert_array_equal(restored_alleles, source)

    expected_phase = raw_block[1]
    if ploidy == 1:
        assert restored_phase in (0, False)
    elif np.all(expected_phase == expected_phase.flat[0]):
        assert bool(restored_phase) == bool(expected_phase.flat[0])
    else:
        np.testing.assert_array_equal(restored_phase, expected_phase)
