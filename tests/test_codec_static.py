import pytest

from gvc import codec


class BrokenPayload:
    def read(self):
        raise AttributeError("synthetic codec read failure")


def test_codec_payload_reader_does_not_swallow_internal_attribute_error():
    with pytest.raises(AttributeError, match="synthetic codec read failure"):
        codec._read_payload(BrokenPayload())


def test_codec_payload_reader_leaves_plain_bytes_unchanged():
    payload = b"abc"
    assert codec._read_payload(payload) is payload
