import argparse

from .binarization import AVAIL_BINARIZATION_MODES
from .codec import AVAIL_CODECS
from .dist import AVAIL_DIST
from .settings import LOG_LEVELS, PROGRAM_DESC, PROGRAM_NAME
from .solver import AVAIL_SOLVERS


def build_parser():
    parser = argparse.ArgumentParser(description=PROGRAM_DESC, prog=PROGRAM_NAME)
    parser.add_argument(
        "-l",
        "--log_level",
        help="log level",
        choices=LOG_LEVELS.keys(),
        default="error",
    )
    subparsers = parser.add_subparsers(dest="mode")

    enc_parser = subparsers.add_parser("encode")
    enc_parser.add_argument("-j", "--num-threads", type=int, default=0)
    enc_parser.add_argument("-b", "--block_size", type=int, default=1024)
    enc_parser.add_argument(
        "--binarization",
        choices=AVAIL_BINARIZATION_MODES,
        default="bit_plane",
    )
    enc_parser.add_argument("--axis", type=int, choices=[0, 1, 2], default=2)
    enc_parser.add_argument("--sort-rows", action="store_true")
    enc_parser.add_argument("--sort-cols", action="store_true")
    enc_parser.add_argument("--dist", choices=AVAIL_DIST, default="ham")
    enc_parser.add_argument("--solver", choices=AVAIL_SOLVERS, default="nn")
    enc_parser.add_argument("--preset-mode", type=int, choices=[0, 1, 2], default=0)
    enc_parser.add_argument("--encoder", choices=AVAIL_CODECS, default=AVAIL_CODECS[0])
    enc_parser.add_argument("input", metavar="input_file_path")
    enc_parser.add_argument("output", metavar="output_file_path")

    dec_parser = subparsers.add_parser("decode")
    dec_parser.add_argument("--pos", type=int, nargs=2)
    dec_parser.add_argument("--samples", type=str)
    dec_parser.add_argument("input", metavar="input_file_path")
    dec_parser.add_argument("output", metavar="output_file_path", nargs="?")

    comp_parser = subparsers.add_parser("compare")
    comp_parser.add_argument("input", metavar="input_file_path")
    comp_parser.add_argument("output", metavar="output_file_path")

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.mode is None:
        parser.print_help()
        return 2

    from . import app

    app.run(args, run_as_module=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
