# import gvc.utils
from .. import utils
from .consts import NROW_LEN, NCOL_LEN

class BinMat():
    def __init__(self, bin_mat):
        self.bin_mat = bin_mat
        
    @staticmethod
    def split_data(data):
        nrows = utils.bytes_to_int(data[:NROW_LEN])
        data = data[NROW_LEN:]
        ncols = utils.bytes_to_int(data[:NCOL_LEN])
        data = data[NCOL_LEN:]
        
        return nrows, ncols, data
    
    def to_bytes(self):

        if self.bin_mat.ndim == 2:
            nrow, ncol = self.bin_mat.shape
        elif self.bin_mat.ndim == 1:
            nrow = self.bin_mat.shape[0]
            ncol = 0

        data_bytes = b''
        data_bytes += utils.int_to_bytes(nrow, NROW_LEN)
        data_bytes += utils.int_to_bytes(ncol, NCOL_LEN)

        bin_mat_bytes = bytes(self.bin_mat.flatten().tolist())

        # data_bytes += gabac.encode(bin_mat_bytes)

        return data_bytes

    @classmethod
    def from_bytes(cls, data):
        raise NotImplementedError()