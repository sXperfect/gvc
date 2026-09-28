import typing as t
BIE_HEADER_LEN = 20

def _get_header(
    bie_header: t.List[int]
) -> t.List[int]:

    ncols = bie_header[4] << 24 | bie_header[5] << 16 | bie_header[6] << 8 | bie_header[7]
    nrows = bie_header[8] << 24 | bie_header[9] << 16 | bie_header[10] << 8 | bie_header[11]

    return nrows, ncols

def get_shape(
    jbig1_bytes:bytes
):

    if isinstance(jbig1_bytes, bytes):
        return _get_header(jbig1_bytes)
    else:
        raise TypeError()