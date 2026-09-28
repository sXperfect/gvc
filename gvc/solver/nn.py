import tspsolve
from .base import *

class NNSolver(AbstractSolver):
    req_dist_mat = True
    
    def __init__(self, preset_mode):
        pass
    
    def sort_rows(self, cost_mat):
        row_ids = tspsolve.nearest_neighbor(cost_mat)
        return row_ids
    
    def sort_cols(self, cost_mat):
        col_ids = tspsolve.nearest_neighbor(cost_mat)
        return col_ids