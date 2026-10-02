import numpy as np
import pytest

from gvc.sort import _sort_matrix


def test_sort_rejects_row_permutation_beyond_uint16_format():
    matrix = np.zeros((np.iinfo(np.uint16).max + 1, 1), dtype=bool)

    with pytest.raises(ValueError, match="row sorting exceeds"):
        _sort_matrix(matrix, sort_row=True)


def test_sort_rejects_column_permutation_beyond_uint16_format():
    matrix = np.zeros((1, np.iinfo(np.uint16).max + 1), dtype=bool)

    with pytest.raises(ValueError, match="column sorting exceeds"):
        _sort_matrix(matrix, sort_col=True)
