from pathlib import Path

import pytest

from gvc.encoder import Encoder


FIXTURE = Path(__file__).parent / "fixtures" / "tiny_diploid.vcf"


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"binarization_name": "missing"}, "unknown binarization"),
        ({"codec_name": "missing"}, "unknown codec"),
        ({"axis": 3}, "axis"),
        ({"block_size": 0}, "block_size"),
        ({"num_threads": -1}, "num_threads"),
        ({"solver": "missing"}, "solver"),
        ({"dist": "missing"}, "distance"),
        ({"preset_mode": 7}, "preset_mode"),
        (
            {"multiprocessing_stall_timeout": 0},
            "multiprocessing_stall_timeout",
        ),
    ],
)
def test_encoder_rejects_invalid_public_arguments(tmp_path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        Encoder(
            FIXTURE,
            tmp_path / "out.gvc",
            **kwargs
        )


def test_encoder_rejects_noninteger_worker_and_block_counts(tmp_path):
    with pytest.raises(TypeError, match="block_size"):
        Encoder(FIXTURE, tmp_path / "out.gvc", block_size=1.5)
    with pytest.raises(TypeError, match="num_threads"):
        Encoder(FIXTURE, tmp_path / "out.gvc", num_threads=True)


def test_encoder_accepts_pathlike_inputs(tmp_path):
    encoder = Encoder(FIXTURE, tmp_path / "out.gvc")
    assert encoder.input_fpath == str(FIXTURE)
    assert encoder.output_fpath == str(tmp_path / "out.gvc")


def test_encoder_rejects_unsupported_input_suffix(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("not a VCF")
    with pytest.raises(ValueError, match=r"\.vcf"):
        Encoder(source, tmp_path / "out.gvc")


def test_encoder_rejects_noncallable_initializer(tmp_path):
    with pytest.raises(TypeError, match="initializer"):
        Encoder(
            FIXTURE,
            tmp_path / "out.gvc",
            multiprocessing_initializer=123,
        )
