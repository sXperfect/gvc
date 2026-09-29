import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


ROOT = Path(__file__).resolve().parents[1]
VCF_FIXTURE = Path(__file__).parent / "fixtures" / "tiny_diploid.vcf"


@pytest.mark.skipif(os.name == "nt", reason="fork-based worker inheritance is POSIX-only")
def test_multiprocessing_encoder_matches_sequential_output(tmp_path):
    script = textwrap.dedent(
        """
        import multiprocessing as mp
        from io import BytesIO
        from pathlib import Path
        import sys

        import numpy as np

        mp.set_start_method("fork", force=True)

        from gvc.codec import MAT_CODECS
        from gvc.codec.jbig import BIE_HEADER_LEN
        from gvc.data_structures.consts import CodecID
        from gvc.decoder import Decoder
        from gvc.encoder import Encoder

        fixture = Path(sys.argv[1])
        root = Path(sys.argv[2])
        sequential = root / "sequential.gvc"
        parallel = root / "parallel.gvc"

        def encode_matrix(matrix):
            matrix = np.asarray(matrix)
            header = bytearray(BIE_HEADER_LEN)
            header[4:8] = int(matrix.shape[1]).to_bytes(4, "big")
            header[8:12] = int(matrix.shape[0]).to_bytes(4, "big")
            payload = BytesIO()
            np.save(payload, matrix, allow_pickle=False)
            return bytes(header) + payload.getvalue()

        def decode_matrix(payload):
            payload = bytes(payload)
            return np.load(BytesIO(payload[BIE_HEADER_LEN:]), allow_pickle=False)

        MAT_CODECS[CodecID.JBIG1]["encoder"] = encode_matrix
        MAT_CODECS[CodecID.JBIG1]["decoder"] = decode_matrix

        common = dict(
            binarization_name="bit_plane",
            axis=2,
            sort_rows=False,
            sort_cols=False,
            transpose=False,
            block_size=1,
            dist="ham",
            solver="nn",
            codec_name="jbig",
            preset_mode=0,
        )

        Encoder(str(fixture), str(sequential), num_threads=0, **common).run()
        Encoder(str(fixture), str(parallel), num_threads=2, **common).run()

        assert sequential.read_bytes() == parallel.read_bytes()

        for name in ("main.npy", "samples.npy", "0.npy", "1.npy", "2.npy"):
            left = Path(str(sequential) + ".metadata") / name
            right = Path(str(parallel) + ".metadata") / name
            np.testing.assert_array_equal(
                np.load(left, allow_pickle=False),
                np.load(right, allow_pickle=False),
            )

        for encoded, name in ((sequential, "seq.txt"), (parallel, "par.txt")):
            decoded = root / name
            decoder = Decoder(str(encoded), str(decoded))
            try:
                decoder.decode()
            finally:
                if decoder._out_f is not None:
                    decoder._out_f.close()
                decoder._f.close()

        expected = "0|1\t1/1\n2/1\t0|2\n./.\t1|0\n"
        assert (root / "seq.txt").read_text() == expected
        assert (root / "par.txt").read_text() == expected
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", script, str(VCF_FIXTURE), str(tmp_path)],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, (
        "multiprocessing regression failed\nstdout:\n{}\nstderr:\n{}".format(
            result.stdout,
            result.stderr,
        )
    )
