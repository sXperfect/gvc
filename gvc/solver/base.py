from abc import ABC, abstractmethod

class AbstractSolver(ABC):
    
    @abstractmethod
    def __init__(self, preset_mode):
        pass
    
    @property
    @abstractmethod
    def req_dist_mat(self):
        pass
    
    @abstractmethod
    def sort_rows(self, mat):
        pass
    
    @abstractmethod
    def sort_cols(self, mat):
        pass