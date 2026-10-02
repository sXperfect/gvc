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
HAPLOID_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_haploid.vcf"
MIXED_PLOIDY_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_mixed_ploidy.vcf"
EXPECTED_GT = "0|1\t1/1\n2/1\t0|2\n./.\t1|0\n"
EXPECTED_HAPLOID_GT = "0\t1\n2\t0\n.\t1\n"
EXPECTED_MIXED_PLOIDY_GT = "0\t1\n2\t0\n0|1\t1/1\n2/1\t0|2\n"


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



def test_encoder_default_codec_name_is_registered(tmp_path):
    encoder = Encoder(
        str(VCF_FIXTURE),
        str(tmp_path / "default.gvc"),
        num_threads=0,
    )
    assert encoder.codec_id == CodecID.JBIG1



def test_complete_haploid_file_roundtrip(framed_test_codec, tmp_path):
    encoded = tmp_path / "haploid.gvc"
    decoded = tmp_path / "haploid.txt"

    Encoder(
        str(HAPLOID_FIXTURE),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=2,
        codec_name="jbig",
        num_threads=0,
    ).run()

    decoder = Decoder(str(encoded), str(decoded))
    try:
        decoder.decode()
    finally:
        if decoder._out_f is not None:
            decoder._out_f.close()
        decoder._f.close()

    assert decoded.read_text() == EXPECTED_HAPLOID_GT



def _close_decoder(decoder):
    if decoder._out_f is not None:
        decoder._out_f.close()
    decoder._f.close()


def test_complete_mixed_ploidy_file_roundtrip(framed_test_codec, tmp_path):
    encoded = tmp_path / "mixed-ploidy.gvc"
    decoded = tmp_path / "mixed-ploidy.txt"

    Encoder(
        str(MIXED_PLOIDY_FIXTURE),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=4,
        codec_name="jbig",
        num_threads=0,
    ).run()

    decoder = Decoder(str(encoded), str(decoded))
    try:
        assert decoder.num_parameter_sets == 2
        assert decoder.num_access_units == 2
        decoder.decode()
    finally:
        _close_decoder(decoder)

    assert decoded.read_text() == EXPECTED_MIXED_PLOIDY_GT


def test_random_access_position_and_sample_subset(framed_test_codec, tmp_path):
    encoded = tmp_path / "random-access.gvc"
    selected = tmp_path / "selected.txt"

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

    decoder = Decoder(str(encoded), str(selected))
    try:
        decoder.random_access([100, 200], "SAMPLE_B")
    finally:
        _close_decoder(decoder)

    assert selected.read_text() == "1/1\n0|2\n"


def test_random_access_sample_only_covers_all_blocks(framed_test_codec, tmp_path):
    encoded = tmp_path / "sample-only.gvc"
    selected = tmp_path / "sample-only.txt"

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

    decoder = Decoder(str(encoded), str(selected))
    try:
        decoder.random_access(None, "SAMPLE_A")
    finally:
        _close_decoder(decoder)

    assert selected.read_text() == "0|1\n2/1\n./.\n"


def test_random_access_empty_interval_writes_nothing(framed_test_codec, tmp_path):
    encoded = tmp_path / "gap.gvc"
    selected = tmp_path / "gap.txt"

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

    decoder = Decoder(str(encoded), str(selected))
    try:
        decoder.random_access([150, 150], None)
    finally:
        _close_decoder(decoder)

    assert selected.read_text() == ""



def _encode_random_access_fixture(framed_test_codec, tmp_path, name):
    encoded = tmp_path / (name + ".gvc")
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
    return encoded


def test_random_access_exact_boundaries_and_multi_block_span(
    framed_test_codec,
    tmp_path,
):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "boundaries"
    )

    cases = [
        ([100, 100], None, "0|1\t1/1\n"),
        ([200, 200], None, "2/1\t0|2\n"),
        ([300, 300], None, "./.\t1|0\n"),
        ([100, 300], None, EXPECTED_GT),
        ([200, 300], "SAMPLE_B", "0|2\n1|0\n"),
    ]
    for index, (position, samples, expected) in enumerate(cases):
        selected = tmp_path / ("selection-{}.txt".format(index))
        decoder = Decoder(str(encoded), str(selected))
        try:
            decoder.random_access(position, samples)
        finally:
            _close_decoder(decoder)
        assert selected.read_text() == expected


def test_random_access_preserves_requested_sample_order(
    framed_test_codec,
    tmp_path,
):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "sample-order"
    )
    selected = tmp_path / "sample-order.txt"

    decoder = Decoder(str(encoded), str(selected))
    try:
        decoder.random_access([100, 100], "SAMPLE_B;SAMPLE_A")
    finally:
        _close_decoder(decoder)

    assert selected.read_text() == "1/1\t0|1\n"


def test_random_access_rejects_invalid_interval_and_unknown_sample(
    framed_test_codec,
    tmp_path,
):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "invalid-random-access"
    )

    decoder = Decoder(str(encoded))
    try:
        with pytest.raises(ValueError, match="start position"):
            decoder.random_access([300, 100], None)
        with pytest.raises(ValueError, match="unknown sample"):
            decoder.random_access([100, 100], "DOES_NOT_EXIST")
    finally:
        _close_decoder(decoder)


def test_random_access_requires_metadata(framed_test_codec, tmp_path):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "missing-metadata"
    )
    metadata = Path(str(encoded) + ".metadata")
    for path in metadata.iterdir():
        path.unlink()
    metadata.rmdir()

    decoder = Decoder(str(encoded))
    try:
        assert decoder.index is None
        with pytest.raises(ValueError, match="metadata"):
            decoder.random_access([100, 100], None)
    finally:
        _close_decoder(decoder)



def test_decoder_rejects_duplicate_parameter_set_id(framed_test_codec, tmp_path):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "duplicate-parameter-set"
    )
    raw = encoded.read_bytes()
    first_len = int.from_bytes(raw[1:5], "big")
    duplicated = raw[:first_len] + raw
    damaged = tmp_path / "duplicate-parameter-set.gvc"
    damaged.write_bytes(duplicated)

    with pytest.raises(ValueError, match="duplicate parameter_set_id"):
        decoder = Decoder(str(damaged))
        _close_decoder(decoder)


def test_decoder_rejects_duplicate_access_unit_id(framed_test_codec, tmp_path):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "duplicate-access-unit-source"
    )
    raw = encoded.read_bytes()
    parameter_len = int.from_bytes(raw[1:5], "big")
    access_start = parameter_len
    access_len = int.from_bytes(raw[access_start + 1:access_start + 5], "big")
    access_unit = raw[access_start:access_start + access_len]
    damaged = tmp_path / "duplicate-access-unit.gvc"
    damaged.write_bytes(
        raw[:access_start + access_len]
        + access_unit
        + raw[access_start + access_len:]
    )

    with pytest.raises(ValueError, match="duplicate access_unit_id"):
        decoder = Decoder(str(damaged))
        _close_decoder(decoder)


def test_decoder_rejects_incomplete_metadata_sidecar(framed_test_codec, tmp_path):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "incomplete-metadata"
    )
    metadata = Path(str(encoded) + ".metadata")
    (metadata / "main.npy").unlink()

    with pytest.raises(FileNotFoundError):
        decoder = Decoder(str(encoded))
        _close_decoder(decoder)



def test_decoder_rejects_metadata_sample_count_mismatch(
    framed_test_codec,
    tmp_path,
):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "sample-count-mismatch"
    )
    metadata = Path(str(encoded) + ".metadata")
    np.save(metadata / "samples.npy", np.array(["ONLY_ONE_SAMPLE"]))

    with pytest.raises(ValueError, match="metadata sample count"):
        decoder = Decoder(str(encoded))
        _close_decoder(decoder)



def test_sequential_encode_failure_preserves_existing_output(
    framed_test_codec,
    monkeypatch,
    tmp_path,
):
    from gvc import encoder as encoder_module

    output = tmp_path / "sequential-existing.gvc"
    output.write_bytes(b"OLD-GVC")
    metadata = Path(str(output) + ".metadata")
    metadata.mkdir()
    (metadata / "marker").write_text("OLD-METADATA")

    def failing_run_core(*args, **kwargs):
        raise RuntimeError("synthetic sequential failure")

    monkeypatch.setattr(encoder_module, "run_core", failing_run_core)

    with pytest.raises(RuntimeError, match="synthetic sequential failure"):
        Encoder(
            str(VCF_FIXTURE),
            str(output),
            binarization_name="bit_plane",
            axis=2,
            sort_rows=False,
            sort_cols=False,
            block_size=2,
            codec_name="jbig",
            num_threads=0,
        ).run()

    assert output.read_bytes() == b"OLD-GVC"
    assert (metadata / "marker").read_text() == "OLD-METADATA"
    assert not list(tmp_path.glob("sequential-existing.gvc.tmp.*"))
    assert not list(tmp_path.glob("sequential-existing.gvc.tmp.*.metadata"))


def test_sequential_success_replaces_existing_output_without_backup_leaks(
    framed_test_codec,
    tmp_path,
):
    output = tmp_path / "sequential-replace.gvc"
    output.write_bytes(b"OLD-GVC")
    metadata = Path(str(output) + ".metadata")
    metadata.mkdir()
    (metadata / "marker").write_text("OLD-METADATA")

    Encoder(
        str(VCF_FIXTURE),
        str(output),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=2,
        codec_name="jbig",
        num_threads=0,
    ).run()

    assert output.read_bytes() != b"OLD-GVC"
    assert not (metadata / "marker").exists()
    assert (metadata / "main.npy").is_file()
    assert (metadata / "samples.npy").is_file()
    assert not list(tmp_path.glob("sequential-replace.gvc.gvc-backup-*"))
    assert not list(tmp_path.glob("sequential-replace.gvc.metadata.gvc-backup-*"))
    assert not list(tmp_path.glob("sequential-replace.gvc.tmp.*"))



def test_decoder_context_manager_closes_owned_files(
    framed_test_codec,
    tmp_path,
):
    encoded = _encode_random_access_fixture(
        framed_test_codec, tmp_path, "decoder-context"
    )
    decoded = tmp_path / "decoder-context.txt"

    with Decoder(str(encoded), str(decoded)) as decoder:
        input_handle = decoder._f
        output_handle = decoder._out_f
        decoder.decode()

    assert input_handle.closed
    assert output_handle.closed
    assert decoder._f is None
    assert decoder._out_f is None



def test_compare_handles_ploidy_induced_block_boundaries(
    framed_test_codec,
    tmp_path,
):
    encoded = tmp_path / "compare-mixed.gvc"

    Encoder(
        str(MIXED_PLOIDY_FIXTURE),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=4,
        codec_name="jbig",
        num_threads=0,
    ).run()

    with Decoder(str(encoded)) as decoder:
        decoder.compare(str(MIXED_PLOIDY_FIXTURE))
