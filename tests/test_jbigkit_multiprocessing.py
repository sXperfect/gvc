import os

import pytest

from gvc.codec import jbigkit
from gvc.decoder import Decoder
from gvc.encoder import Encoder


def _require_real_jbig():
    if os.environ.get("GVC_REQUIRE_JBIGKIT") != "1":
        pytest.skip("real JBIG-KIT multiprocessing stress is release-gate only")
    if not jbigkit.available():
        pytest.fail("real JBIG-KIT release dependency is unavailable")


def _write_stress_vcf(path, variants=16):
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    ]
    patterns = [
        ("0|1", "1/1", "0/0", "1|0"),
        ("2/1", "0|2", "1/2", "0/0"),
        ("./.", "1|0", "0/1", "1/1"),
        ("1|1", "0/0", "0|1", "1/0"),
    ]
    lines = list(header)
    expected = []
    for index in range(variants):
        genotypes = patterns[index % len(patterns)]
        pos = 100 + index * 10
        lines.append(
            "1\t{}\t.\tA\tC,G\t60\tPASS\t.\tGT\t{}".format(
                pos,
                "\t".join(genotypes),
            )
        )
        expected.append("\t".join(genotypes))
    path.write_text("\n".join(lines) + "\n")
    return "\n".join(expected) + "\n"


def _close_decoder(decoder):
    if decoder._out_f is not None:
        decoder._out_f.close()
    decoder._f.close()


def test_real_jbigkit_four_worker_queue_pressure(tmp_path):
    _require_real_jbig()
    jbigkit.install()

    source = tmp_path / "stress.vcf"
    expected = _write_stress_vcf(source, variants=16)
    encoded = tmp_path / "stress.gvc"
    decoded = tmp_path / "decoded.txt"

    Encoder(
        str(source),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=1,
        codec_name="jbig",
        preset_mode=0,
        num_threads=4,
        multiprocessing_stall_timeout=30,
        multiprocessing_initializer=jbigkit.install,
    ).run()

    decoder = Decoder(str(encoded), str(decoded))
    try:
        decoder.decode()
    finally:
        _close_decoder(decoder)

    assert decoded.read_text() == expected
    metadata = tmp_path / "stress.gvc.metadata"
    assert (metadata / "main.npy").is_file()
    assert not list(tmp_path.glob("stress.gvc.tmp.*"))


def test_real_jbigkit_parallel_random_access_after_stress(tmp_path):
    _require_real_jbig()
    jbigkit.install()

    source = tmp_path / "stress-ra.vcf"
    _write_stress_vcf(source, variants=8)
    encoded = tmp_path / "stress-ra.gvc"
    selected = tmp_path / "selected.txt"

    Encoder(
        str(source),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=1,
        codec_name="jbig",
        preset_mode=0,
        num_threads=2,
        multiprocessing_stall_timeout=30,
        multiprocessing_initializer=jbigkit.install,
    ).run()

    decoder = Decoder(str(encoded), str(selected))
    try:
        decoder.random_access([120, 150], "S4;S2")
    finally:
        _close_decoder(decoder)

    assert selected.read_text() == "1/1\t1|0\n1/0\t0/0\n1|0\t1/1\n0/0\t0|2\n"
