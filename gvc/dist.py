import numpy as np
from scipy.spatial.distance import pdist, squareform
from numba import jit
# from numba import jit, uint16, boolean

@jit(nopython=True, cache=True)
#@jit(uint16(boolean[:], boolean[:]), nopython=True, cache=True)
def hamming_rl_dist(vector1, vector2):
    """
    Compute Hamming - Run Length distance given 2 vectors.

    Parameters
    ----------
    vector1 : ndarray (1d)
        First vector, contains binary values
    vector1 : ndarray (1d)
        Second vector, contains binary values

    Returns
    -------
    cost : int
        Cost given 2 vector using hamming run length distance

    TODO: Fix example
    For example:
    0 1 2 3 4 5 6 7\n 
    0|0|1|0|1|0|0|0\n
    0|0|0|1|0|1|0|0\n
        d d d d    \n
    diff_locs : 2,3,4,5\n
    diff = 1 (entry)\n

    0 1 2 3 4 5 6 7\n
    0|1|0|0|0|1|0|0\n
    0|0|1|0|0|0|1|0\n
      d d     d d  \n
    diff_locs : 1,2,5,6\n
    diff : 2 (entries)\n

    """
    diff_locs = np.where(np.logical_not(vector1 == vector2))[0]

    if np.any(diff_locs):
        return np.sum(np.diff(diff_locs) > 1) + 1
    else:
        return 0

DIST_FUNC = {
    "ham": "hamming",
    "ham_rl": hamming_rl_dist
}

AVAIL_DIST = [k for k in DIST_FUNC.keys()]

def comp_cost_mat(
    bin_mat, 
    dist_f:str
):
    """_summary_

    Args:
        bin_mat (array_like): An m by n array of m original observations in an n-dimensional space
        dist_f (str): name of distance function

    Returns:
        np.array: distance matrix
    """
    
    dist_mat = squareform(
        pdist(bin_mat, metric=DIST_FUNC[dist_f])
    )
    
    return dist_mat