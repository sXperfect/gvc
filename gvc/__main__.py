import argparse

import gvc.encoder
import gvc.decoder
import gvc.settings
import gvc.binarization
import gvc.sort
import gvc.app

from .settings import LOG_LEVELS
from .dist import AVAIL_DIST
from .solver import AVAIL_SOLVERS
from .codec import AVAIL_CODECS

block_size_help = 'block size (in number of records), will be ignored if -d/--decode is used'

parser = argparse.ArgumentParser(description=gvc.settings.PROGRAM_DESC, prog=gvc.settings.PROGRAM_NAME)
parser.add_argument('-l', '--log_level', help='log level', choices=LOG_LEVELS.keys(), default='info')
parser.add_argument('-c', '--compare', help='compare', action='store_true')
parser.add_argument('-j', '--num-threads', help='?', type=int, default=0)
parser.add_argument('-b', '--block_size', help=block_size_help, type=int, default=1024)
parser.add_argument('--binarization', help='binarization method', choices=gvc.binarization.binarization_str_to_flag.keys(), default='bit_plane')
parser.add_argument('--axis', help='The axis along which the matrices of bit_plane will be joined.', type=int, choices=[0, 1, 2], default=2)
parser.add_argument('--sort-rows', help='Sort rows', action='store_true')
parser.add_argument('--sort-cols', help='Sort columns', action='store_true')
parser.add_argument('-t', '--transpose', help='Tranpose binary matrix', action='store_true')
parser.add_argument('--dist', help='Distance to compute cost matrix', choices=AVAIL_DIST, default='ham')
parser.add_argument('--solver', help='Traveling Salesman Problem solver', choices=AVAIL_SOLVERS, default='nn')
parser.add_argument('--preset-mode', type=int, choices=[0, 1, 2], default=0)
parser.add_argument('--encoder', help='?', choices=AVAIL_CODECS, default='jbig')
parser.add_argument('mode', choices=["encode", "decode", "stat", "compare", "cprofile"])
parser.add_argument('input', metavar='input_file_path', help='input file path')
parser.add_argument('output', metavar='output_file_path', help='output file path', nargs='?')
args = parser.parse_args()

gvc.app.run(args, run_as_module=False)
