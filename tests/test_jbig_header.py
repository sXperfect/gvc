from gvc.codec import jbig


def test_jbig_shape_header_validation():
    import pytest

    with pytest.raises(ValueError, match="truncated"):
        jbig.get_shape(b"\x00" * (jbig.BIE_HEADER_LEN - 1))

    header = bytearray(jbig.BIE_HEADER_LEN)
    with pytest.raises(ValueError, match="non-positive"):
        jbig.get_shape(bytes(header))

    header[4:8] = (7).to_bytes(4, "big")
    header[8:12] = (5).to_bytes(4, "big")
    assert jbig.get_shape(bytes(header)) == (5, 7)
