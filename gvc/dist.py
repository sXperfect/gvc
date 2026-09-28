"""Distance functions used by the sorting stage."""

import numpy as np
from scipy.spatial.distance import pdist, squareform

try:
    from numba import jit
except ImportError:
    def jit(*args, **kwargs):
        if args and callable(args[0]) and len(args) == 1 and not kwargs:
            return args[0]

        def decorate(func):
            return func

        return decorate


@jit(nopython=True, cache=True)
def hamming_rl_dist(vector1, vector2):
    """Return the number of contiguous mismatch runs between two vectors."""
    diff_locs = np.where(np.logical_not(vector1 == vector2))[0]
    if np.any(diff_locs):
        return np.sum(np.diff(diff_locs) > 1) + 1
    return 0


DIST_FUNC = {
    "ham": "hamming",
    "ham_rl": hamming_rl_dist,
}

AVAIL_DIST = list(DIST_FUNC)


def comp_cost_mat(bin_mat, dist_f):
    """Compute a square pairwise cost matrix for rows of bin_mat."""
    if dist_f not in DIST_FUNC:
        raise ValueError("unknown distance function: {}".format(dist_f))
    return squareform(pdist(bin_mat, metric=DIST_FUNC[dist_f]))
