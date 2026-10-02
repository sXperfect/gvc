import os

import numpy as np
import pytest

from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.decoder import decode_encoded_variants
from gvc.encoder import run_core
from gvc.reader import vcf_genotypes_reader


FIXTURE_ENV = "GVC_HISTORICAL_FIXTURE"
MAX_BLOCKS_ENV = "GVC_HISTORICAL_MAX_BLOCKS"


def _historical_fixture():
    path = os.environ.get(FIXTURE_ENV)
    if not path:
        pytest.skip("historical LUH fixture is release-gate only")
    if not os.path.isfile(path):
        pytest.fail("{} does not exist: {}".format(FIXTURE_ENV, path))
    return path


@pytest.mark.parametrize(
    "binarization_id,axis",
    [
        (BinarizationID.BIT_PLANE, 2),
        (BinarizationID.ROW_BIN_SPLIT, 0),
    ],
)
def test_original_luh_vcf_fixture_core_roundtrip(
    binarization_id,
    axis,
):
    fixture = _historical_fixture()
    max_blocks = int(os.environ.get(MAX_BLOCKS_ENV, "1"))
    if max_blocks <= 0:
        raise ValueError("{} must be positive".format(MAX_BLOCKS_ENV))

    ps_params = [
        binarization_id,
        CodecID.JBIG1,
        axis,
        False,
        False,
        False,
    ]
    tsp_params = ["ham", "nn", 0]

    seen = 0
    for raw_block in vcf_genotypes_reader(fixture, None, block_size=2048):
        original_alleles = raw_block[0].copy()
        original_phase = raw_block[1].copy()

        block, parameter_set = run_core(raw_block, ps_params, tsp_params)
        restored_alleles, restored_phase = decode_encoded_variants(
            parameter_set,
            block.block_payload,
            ret_gt=False,
        )

        np.testing.assert_array_equal(restored_alleles, original_alleles)

        if parameter_set.p == 1:
            assert restored_phase in (0, False)
        elif parameter_set.encode_phase_data:
            np.testing.assert_array_equal(restored_phase, original_phase)
        elif original_phase.size:
            assert bool(restored_phase) == bool(original_phase.flat[0])

        seen += 1
        if seen >= max_blocks:
            break

    assert seen == max_blocks
