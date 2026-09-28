"""Deterministic nearest-neighbor ordering used by GVC sorting."""

import numpy as np

from .base import AbstractSolver


def nearest_neighbor(cost_mat):
    """Return a deterministic nearest-neighbor tour starting at index 0.

    Ties are resolved by the lowest unvisited index. The function intentionally
    implements only the small contract GVC needs and replaces the historical
    dependency on the unmaintained tspsolve package.
    """
    costs = np.asarray(cost_mat)
    if costs.ndim != 2 or costs.shape[0] != costs.shape[1]:
        raise ValueError("cost matrix must be square")

    n_items = costs.shape[0]
    if n_items == 0:
        return np.empty(0, dtype=np.int64)

    route = np.empty(n_items, dtype=np.int64)
    unvisited = np.ones(n_items, dtype=bool)
    current = 0

    for pos in range(n_items):
        route[pos] = current
        unvisited[current] = False
        if pos == n_items - 1:
            break

        candidates = np.flatnonzero(unvisited)
        next_offset = int(np.argmin(costs[current, candidates]))
        current = int(candidates[next_offset])

    return route


class NNSolver(AbstractSolver):
    req_dist_mat = True

    def __init__(self, preset_mode):
        self.preset_mode = preset_mode

    def sort_rows(self, cost_mat):
        return nearest_neighbor(cost_mat)

    def sort_cols(self, cost_mat):
        return nearest_neighbor(cost_mat)
