from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from gvc.codec import MAT_CODECS
from gvc.codec.jbig import BIE_HEADER_LEN
from gvc.data_structures.consts import CodecID
from gvc.decoder import Decoder
from gvc.encoder import Encoder


VCF_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_diploid.vcf"
EXPECTED_GT = "0|1\t1/1\n2/1\t0|2\n./.\t1|0\n"


def _framed_array_encode(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError("test codec only accepts matrices")

    header = bytearray(BIE_HEADER_LEN)
    header[4:8] = int(matrix.shape[1]).to_bytes(4, "big")
    header[8:12] = int(matrix.shape[0]).to_bytes(4, "big")

    payload = BytesIO()
    np.save(payload, matrix, allow_pickle=False)
    return bytes(header) + payload.getvalue()


def _framed_array_decode(payload):
    payload = bytes(payload)
    if len(payload) < BIE_HEADER_LEN:
        raise ValueError("test codec payload is truncated")
    return np.load(BytesIO(payload[BIE_HEADER_LEN:]), allow_pickle=False)


@pytest.fixture
def framed_test_codec(monkeypatch):
    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", _framed_array_encode)
    monkeypatch.setitem(codec, "decoder", _framed_array_decode)


@pytest.mark.parametrize(
    "binarization,axis",
    [
        ("bit_plane", 2),
        ("row_bin_split", 0),
    ],
)
def test_complete_file_encode_decode_roundtrip(
    framed_test_codec,
    tmp_path,
    binarization,
    axis,
):
    encoded = tmp_path / ("roundtrip-{}.gvc".format(binarization))
    decoded = tmp_path / ("roundtrip-{}.txt".format(binarization))

    Encoder(
        str(VCF_FIXTURE),
        str(encoded),
        binarization_name=binarization,
        axis=axis,
        sort_rows=False,
        sort_cols=False,
        transpose=False,
        block_size=2,
        dist="ham",
        solver="nn",
        codec_name="jbig",
        preset_mode=0,
        num_threads=0,
    ).run()

    assert encoded.is_file()
    assert encoded.stat().st_size > 0
    metadata = Path(str(encoded) + ".metadata")
    assert (metadata / "main.npy").is_file()
    assert (metadata / "samples.npy").is_file()

    decoder = Decoder(str(encoded), str(decoded))
    try:
        assert decoder.num_parameter_sets >= 1
        assert decoder.num_access_units >= 1
        decoder.decode()
    finally:
        if decoder._out_f is not None:
            decoder._out_f.close()
        decoder._f.close()

    assert decoded.read_text() == EXPECTED_GT


def test_decoder_rejects_truncated_complete_file(framed_test_codec, tmp_path):
    encoded = tmp_path / "full.gvc"
    Encoder(
        str(VCF_FIXTURE),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=2,
        codec_name="jbig",
        num_threads=0,
    ).run()

    raw = encoded.read_bytes()
    for cut in (1, 4, len(raw) // 2, len(raw) - 1):
        damaged = tmp_path / ("truncated-{}.gvc".format(cut))
        damaged.write_bytes(raw[:cut])
        with pytest.raises((EOFError, ValueError, TypeError)):
            decoder = Decoder(str(damaged))
            decoder._f.close()
