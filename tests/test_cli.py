import subprocess
import sys

import pytest

from gvc import app

from gvc.__main__ import build_parser, main


def test_encode_parser_defaults():
    args = build_parser().parse_args(["encode", "input.vcf", "output.gvc"])
    assert args.mode == "encode"
    assert args.block_size == 1024
    assert args.dist == "ham"
    assert args.solver == "nn"


def test_no_mode_returns_usage_error(capsys):
    assert main([]) == 2
    captured = capsys.readouterr()
    assert "usage:" in captured.out



def test_module_cli_help_runs_in_subprocess():
    result = subprocess.run(
        [sys.executable, "-m", "gvc", "--help"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    assert "usage:" in result.stdout.lower()



def test_cli_decode_failure_preserves_existing_output(monkeypatch, tmp_path):
    output = tmp_path / "decoded.txt"
    output.write_text("OLD")

    class FailingDecoder:
        def __init__(self, input_fpath, output_fpath=None):
            self.output_fpath = output_fpath

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def decode(self):
            with open(self.output_fpath, "w") as handle:
                handle.write("PARTIAL")
            raise RuntimeError("synthetic decode failure")

    monkeypatch.setattr(app, "Decoder", FailingDecoder)

    with pytest.raises(RuntimeError, match="synthetic decode failure"):
        app._run_decode_transaction(
            "input.gvc",
            str(output),
            lambda decoder: decoder.decode(),
        )

    assert output.read_text() == "OLD"
    assert not list(tmp_path.glob("decoded.txt.tmp.*"))


def test_cli_decode_success_atomically_replaces_output(monkeypatch, tmp_path):
    output = tmp_path / "decoded.txt"
    output.write_text("OLD")

    class SuccessfulDecoder:
        def __init__(self, input_fpath, output_fpath=None):
            self.output_fpath = output_fpath

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def decode(self):
            with open(self.output_fpath, "w") as handle:
                handle.write("NEW")

    monkeypatch.setattr(app, "Decoder", SuccessfulDecoder)

    app._run_decode_transaction(
        "input.gvc",
        str(output),
        lambda decoder: decoder.decode(),
    )

    assert output.read_text() == "NEW"
    assert not list(tmp_path.glob("decoded.txt.tmp.*"))
