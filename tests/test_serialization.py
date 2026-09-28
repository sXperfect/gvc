from io import BytesIO

import numpy as np

from gvc.bitstream import BitstreamReader
from gvc.codec import MAT_CODECS
from gvc.common import create_parameter_set
from gvc.data_structures import AccessUnit, Block, ParameterSet
from gvc.data_structures.consts import BinarizationID, CodecID, DataUnitType
from gvc.data_structures.data_unit import DataUnitHeader
from gvc.decoder import decode_encoded_variants
from gvc.encoder import run_core


def _array_encode(matrix):
    buffer = BytesIO()
    np.save(buffer, np.asarray(matrix), allow_pickle=False)
    return buffer.getvalue()


def _array_decode(payload):
    return np.load(BytesIO(payload), allow_pickle=False)


def _install_test_codec(monkeypatch):
    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", _array_encode)
    monkeypatch.setitem(codec, "decoder", _array_decode)


def test_parameter_set_binary_roundtrip():
    phase = np.array([[False, True], [True, False]], dtype=bool)
    parameter_set = create_parameter_set(
        missing_rep_val=4,
        na_rep_val=None,
        p=2,
        phasing_matrix=phase,
        additional_info=2,
        binarization_id=BinarizationID.BIT_PLANE,
        codec_id=CodecID.JBIG1,
        axis=2,
        sort_rows=True,
        sort_cols=False,
    )

    raw = parameter_set.to_bytes()
    reader = BitstreamReader(BytesIO(raw))
    data_unit_type = reader.read_bytes(1, ret_int=True)
    assert data_unit_type == DataUnitType.PARAMETER_SET

    header = DataUnitHeader.from_bitstream(data_unit_type, reader)
    restored = ParameterSet.from_bitstream(reader, header)

    assert restored == parameter_set


def test_block_and_access_unit_binary_roundtrip(monkeypatch):
    _install_test_codec(monkeypatch)

    allele_matrix = np.array(
        [
            [0, 1, 2, 0],
            [2, 1, 0, 3],
            [1, 0, 3, 2],
        ],
        dtype=np.uint8,
    )
    phase_matrix = np.array(
        [
            [False, True],
            [True, False],
            [False, False],
        ],
        dtype=bool,
    )

    block, parameter_set = run_core(
        (allele_matrix, phase_matrix, 2, None, None),
        [
            BinarizationID.BIT_PLANE,
            CodecID.JBIG1,
            2,
            True,
            True,
            False,
        ],
        ["ham", "nn", 0],
    )

    block_reader = BitstreamReader(BytesIO(block.to_bytes()))
    restored_block = Block.from_bitstream(block_reader, parameter_set)
    restored_alleles, restored_phase = decode_encoded_variants(
        parameter_set,
        restored_block.block_payload,
        ret_gt=False,
    )
    np.testing.assert_array_equal(restored_alleles, allele_matrix)
    np.testing.assert_array_equal(restored_phase, phase_matrix)

    access_unit = AccessUnit.from_blocks(0, parameter_set.parameter_set_id, [block])
    access_reader = BitstreamReader(BytesIO(access_unit.to_bytes()))
    data_unit_type = access_reader.read_bytes(1, ret_int=True)
    assert data_unit_type == DataUnitType.ACCESS_UNIT

    restored_access = AccessUnit.from_bitstream(
        access_reader,
        {parameter_set.parameter_set_id: parameter_set},
    )
    assert restored_access.header.access_unit_id == 0
    assert restored_access.num_blocks == 1

    restored_alleles, restored_phase = decode_encoded_variants(
        parameter_set,
        restored_access.blocks[0].block_payload,
        ret_gt=False,
    )
    np.testing.assert_array_equal(restored_alleles, allele_matrix)
    np.testing.assert_array_equal(restored_phase, phase_matrix)
