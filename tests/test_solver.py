import numpy as np
import pytest

from gvc.solver.nn import nearest_neighbor


def test_nearest_neighbor_follows_lowest_cost_path():
    costs = np.array(
        [
            [0.0, 1.0, 4.0],
            [1.0, 0.0, 2.0],
            [4.0, 2.0, 0.0],
        ]
    )
    np.testing.assert_array_equal(nearest_neighbor(costs), [0, 1, 2])


def test_nearest_neighbor_ties_are_deterministic():
    costs = np.ones((4, 4), dtype=float)
    np.fill_diagonal(costs, 0.0)
    np.testing.assert_array_equal(nearest_neighbor(costs), [0, 1, 2, 3])


def test_nearest_neighbor_accepts_empty_matrix():
    route = nearest_neighbor(np.empty((0, 0)))
    assert route.shape == (0,)


def test_nearest_neighbor_rejects_non_square_matrix():
    with pytest.raises(ValueError, match="square"):
        nearest_neighbor(np.zeros((2, 3)))
