from types import SimpleNamespace

import numpy as np
import pytest

from gvc.data_structures.bin_mat import BinMat
from gvc.data_structures.index import Index


def _decoder_context(blocks, parameter_set_id=4):
    access_unit = SimpleNamespace(
        blocks=list(blocks),
        header=SimpleNamespace(parameter_set_id=parameter_set_id),
    )
    return SimpleNamespace(access_units={0: access_unit})


@pytest.mark.parametrize(
    "matrix",
    [
        np.array([True, False, True], dtype=bool),
        np.array([[0, 1, 0], [1, 1, 0]], dtype=np.uint8),
        np.zeros((0, 3), dtype=np.uint8),
    ],
)
def test_binmat_roundtrip(matrix):
    encoded = BinMat(matrix).to_bytes()
    restored = BinMat.from_bytes(encoded).bin_mat
    np.testing.assert_array_equal(restored, matrix.astype(bool))


def test_binmat_rejects_nonbinary_and_truncated_payloads():
    with pytest.raises(ValueError, match="0 or 1"):
        BinMat(np.array([0, 2], dtype=np.uint8))

    valid = BinMat(np.array([[0, 1], [1, 0]], dtype=np.uint8)).to_bytes()
    with pytest.raises(ValueError, match="length mismatch"):
        BinMat.from_bytes(valid[:-1])


def test_index_overlap_queries_and_sample_lookup(tmp_path):
    np.save(
        tmp_path / "main.npy",
        np.array([[100, 199], [300, 399], [500, 599]], dtype=np.uint64),
    )
    np.save(tmp_path / "samples.npy", np.array(["S1", "S2", "S3"]))
    blocks = [object(), object(), object()]
    index = Index(str(tmp_path), _decoder_context(blocks))

    assert index.query_blk(1, 50).shape == (0, 3)
    assert index.query_blk(150, 350)[:, 0].tolist() == [0, 1]
    assert index.query_blk(250, 275).shape == (0, 3)
    assert index.query_blk(600, 700).shape == (0, 3)

    np.testing.assert_array_equal(
        index.query_columns("S3;S1"),
        np.array([2, 0], dtype=np.uint32),
    )
    with pytest.raises(ValueError, match="unknown sample"):
        index.query_columns("S4")


def test_index_validates_metadata_and_row_positions(tmp_path):
    np.save(tmp_path / "main.npy", np.array([[100, 120]], dtype=np.uint64))
    np.save(tmp_path / "samples.npy", np.array(["S1"]))
    np.save(tmp_path / "0.npy", np.array([100, 110, 120], dtype=np.uint64))

    index = Index(str(tmp_path), _decoder_context([object()]))
    row_slice = index.get_row_mask(0, 105, 120)
    assert (row_slice.start, row_slice.stop) == (1, 3)

    with pytest.raises(ValueError, match="start position"):
        index.query_blk(20, 10)
    with pytest.raises(IndexError, match="block_id"):
        index.get_row_mask(1, 100, 110)


def test_index_rejects_metadata_block_count_mismatch(tmp_path):
    np.save(
        tmp_path / "main.npy",
        np.array([[100, 199], [300, 399]], dtype=np.uint64),
    )
    np.save(tmp_path / "samples.npy", np.array(["S1"]))
    with pytest.raises(ValueError, match="block count"):
        Index(str(tmp_path), _decoder_context([object()]))



def test_index_rejects_noninteger_root_positions(tmp_path):
    np.save(tmp_path / "main.npy", np.array([[100.0, 120.0]], dtype=np.float64))
    np.save(tmp_path / "samples.npy", np.array(["S1"]))

    with pytest.raises(TypeError, match="root index positions"):
        Index(str(tmp_path), _decoder_context([object()]))


def test_index_rejects_block_positions_that_disagree_with_root(tmp_path):
    np.save(tmp_path / "main.npy", np.array([[100, 120]], dtype=np.uint64))
    np.save(tmp_path / "samples.npy", np.array(["S1"]))
    np.save(tmp_path / "0.npy", np.array([100, 110, 119], dtype=np.uint64))

    index = Index(str(tmp_path), _decoder_context([object()]))
    with pytest.raises(ValueError, match="does not match root index bounds"):
        index.get_row_mask(0, 100, 120)
