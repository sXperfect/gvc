from queue import Queue

import pytest

from gvc.codec import MAT_CODECS
from gvc.data_structures.consts import CodecID
from gvc.decoder import Decoder
from gvc.encoder import Encoder, worker_writer

from tests.mp_test_codec import decode as mp_decode
from tests.mp_test_codec import encode as mp_encode
from tests.test_file_roundtrip import (
    EXPECTED_GT,
    VCF_FIXTURE,
    _close_decoder,
)


def test_multiprocessing_encoder_matches_single_process(monkeypatch, tmp_path):
    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", mp_encode)
    monkeypatch.setitem(codec, "decoder", mp_decode)

    single = tmp_path / "single.gvc"
    multi = tmp_path / "multi.gvc"

    common = dict(
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        transpose=False,
        block_size=1,
        codec_name="jbig",
        preset_mode=0,
    )
    Encoder(str(VCF_FIXTURE), str(single), num_threads=0, **common).run()
    Encoder(str(VCF_FIXTURE), str(multi), num_threads=2, **common).run()

    # Worker scheduling must not affect serialized block order or metadata.
    assert multi.read_bytes() == single.read_bytes()

    decoder = Decoder(str(multi), str(tmp_path / "decoded.txt"))
    try:
        decoder.decode()
    finally:
        _close_decoder(decoder)
    assert (tmp_path / "decoded.txt").read_text() == EXPECTED_GT


def test_writer_orders_out_of_order_worker_results(monkeypatch, tmp_path):
    from gvc import encoder as encoder_module

    queue = Queue()
    stored = []

    class FakeParameterSet:
        parameter_set_id = 0

        def __eq__(self, other):
            return isinstance(other, FakeParameterSet)

        def to_bytes(self):
            return b"P"

    class FakeBlock:
        def __init__(self, value):
            self.value = value

        def __len__(self):
            return 1

    parameter_set = FakeParameterSet()

    def fake_store_access_unit(output_f, access_unit_id, param_set, blocks):
        values = [block.value for block in blocks]
        stored.extend(values)
        output_f.write(bytes(values))

    monkeypatch.setattr(
        encoder_module.gvc.common,
        "store_access_unit",
        fake_store_access_unit,
    )

    # Deliberately emulate workers finishing 2, 0, 1.
    queue.put((2, parameter_set, FakeBlock(2)))
    queue.put((0, parameter_set, FakeBlock(0)))
    queue.put((1, parameter_set, FakeBlock(1)))
    queue.put(None)
    queue.put(None)

    output = tmp_path / "ordered.gvc"
    worker_writer(queue, str(output), num_processes=2)

    assert stored == [0, 1, 2]
    assert output.read_bytes() == b"P\x00\x01\x02"


def test_writer_rejects_missing_block(tmp_path):
    queue = Queue()
    parameter_set = type(
        "FakeParameterSet",
        (),
        {
            "parameter_set_id": 0,
            "to_bytes": lambda self: b"P",
            "__eq__": lambda self, other: True,
        },
    )()
    block = type("FakeBlock", (), {"__len__": lambda self: 1})()

    queue.put((1, parameter_set, block))
    queue.put(None)

    with pytest.raises(RuntimeError, match="missing encoded block"):
        worker_writer(queue, str(tmp_path / "missing.gvc"), num_processes=1)
