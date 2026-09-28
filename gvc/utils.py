import os
from time import perf_counter


def line_cnt(fpath):
    with open(fpath) as f:
        return sum(1 for _ in f)


def int2bstr(val, len_in_byte, order="big"):
    return int(val).to_bytes(len_in_byte, order)


def bstr2int(data, order="big"):
    return int.from_bytes(data, order)


def int_to_bytes(val, len_in_byte, order="big"):
    return int2bstr(val, len_in_byte, order)


def bytes_to_int(data, order="big"):
    return bstr2int(data, order)


def check_executable(path):
    if not os.path.isfile(path):
        raise FileNotFoundError("this is not a file: {}".format(path))
    if not os.access(path, os.X_OK):
        raise FileNotFoundError("file is not executable: {}".format(path))


class catchtime:
    def __enter__(self):
        self.time = perf_counter()
        return self

    def __exit__(self, exc_type, value, traceback):
        self.time = perf_counter() - self.time
