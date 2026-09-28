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
