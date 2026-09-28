import os
import logging as log
from time import perf_counter

def line_cnt(fpath):
    with open(fpath) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# def int_to_bytes(val, len_in_byte, order='big'):
#     assert isinstance(val, int)
#     return (val).to_bytes(len_in_byte, order)

# def bytes_to_int(data, order='big'):
#     return int.from_bytes(data, order)

def int2bstr(val, len_in_byte, order='big'):
    return int(val).to_bytes(len_in_byte, order) 

def bstr2int(data, order='big'):
    return int.from_bytes(data, order)

def check_executable(path):
    if not os.path.isfile(path):
        log.error("this is not a file: {}".format(path))
        raise FileNotFoundError("this is not a file: {}".format(path))
    if not os.access(path, os.X_OK):
        raise FileNotFoundError("file is not executable: {}".format(path))

class catchtime:
    def __enter__(self):
        self.time = perf_counter()
        return self

    def __exit__(self, type, value, traceback):
        self.time = perf_counter() - self.time
        # self.readout = f'Time: {self.time:.3f} seconds'