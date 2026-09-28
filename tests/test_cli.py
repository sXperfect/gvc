import subprocess
import sys

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
