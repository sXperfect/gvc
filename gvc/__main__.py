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
from .binarization import AVAIL_BINARIZATION_MODES

block_size_help = 'block size (in number of records), will be ignored if -d/--decode is used'

parser = argparse.ArgumentParser(description=gvc.settings.PROGRAM_DESC, prog=gvc.settings.PROGRAM_NAME)
parser.add_argument('-l', '--log_level', help='log level', choices=LOG_LEVELS.keys(), default='info')
subparsers = parser.add_subparsers(dest='mode')

#? Encode
encode_parser = subparsers.add_parser('encode')
encode_parser.add_argument('-l', '--log_level', help='log level', choices=LOG_LEVELS.keys(), default='info')
encode_parser.add_argument('-c', '--compare', help='compare', action='store_true')
encode_parser.add_argument('-j', '--num-threads', help='?', type=int, default=0)
encode_parser.add_argument('-b', '--block_size', help=block_size_help, type=int, default=1024)
encode_parser.add_argument('--binarization', help='binarization method', choices=AVAIL_BINARIZATION_MODES, default='bit_plane')
encode_parser.add_argument('--axis', help='The axis along which the matrices of bit_plane will be joined.', type=int, choices=[0, 1, 2], default=2)
encode_parser.add_argument('--sort-rows', help='Sort rows', action='store_true')
encode_parser.add_argument('--sort-cols', help='Sort columns', action='store_true')
encode_parser.add_argument('-t', '--transpose', help='Tranpose binary matrix', action='store_true')
encode_parser.add_argument('--dist', help='Distance to compute cost matrix', choices=AVAIL_DIST, default='ham')
encode_parser.add_argument('--solver', help='Traveling Salesman Problem solver', choices=AVAIL_SOLVERS, default='nn')
encode_parser.add_argument('--preset-mode', type=int, choices=[0, 1, 2], default=0)
encode_parser.add_argument('--encoder', help='?', choices=AVAIL_CODECS, default='jbig')
encode_parser.add_argument('input', metavar='input_file_path', help='input file path')
encode_parser.add_argument('output', metavar='output_file_path', help='output file path')

#? Decode
decode_parser = subparsers.add_parser('decode')
decode_parser.add_argument('input', metavar='input_file_path', help='input file path')
decode_parser.add_argument('output', metavar='output_file_path', help='output file path')

#? Random Access
ra_parser = subparsers.add_parser('random-access')
ra_parser.add_argument('input', metavar='input_file_path', help='input file path')
ra_parser.add_argument('pos', type=int, help='Position of genomic variants to be decoded. If "--end END_POS" is set, then all variants at given range are decoded')
ra_parser.add_argument('--end', type=int, default=-1, help='End position')
ra_parser.add_argument('--samples', default=None, help='List of samples to be decoded. Default: all samples')

args = parser.parse_args()

gvc.app.run(args, run_as_module=False)
