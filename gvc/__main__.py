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

num_threads_help = 'The number of threads'
block_size_help = 'Block size (in number of records), will be ignored if -d/--decode is used'
binarization_help = 'Binarization method'
axis_help = 'Axis, along which the binary matrices of bit_plane will be concatenated. 0 for first dimension (row), 1 for second dimension (column) and 2 for no concatenation.'
preset_mode_help = 'Preset mode of the traveling salesman problem solver'
encoder_help = 'Entropy codec'

parser = argparse.ArgumentParser(description=gvc.settings.PROGRAM_DESC, prog=gvc.settings.PROGRAM_NAME)
parser.add_argument('-l', '--log_level', help='log level', choices=LOG_LEVELS.keys(), default='info')
subparsers = parser.add_subparsers(dest='mode')
# parser.add_argument('input', metavar='input_file_path', help='input file path')
# parser.add_argument('output', metavar='output_file_path', help='output file path', nargs='?')

#? Encode
enc_parser = subparsers.add_parser('encode')
enc_parser.add_argument('-l', '--log_level', help='log level', choices=LOG_LEVELS.keys(), default='info')
# enc_parser.add_argument('-c', '--compare', help='compare', action='store_true')
enc_parser.add_argument('-j', '--num-threads', help=num_threads_help, type=int, default=0)
enc_parser.add_argument('-b', '--block_size', help=block_size_help, type=int, default=1024)
enc_parser.add_argument('--binarization', help=binarization_help, choices=AVAIL_BINARIZATION_MODES, default='bit_plane')
enc_parser.add_argument('--axis', type=int, choices=[0, 1, 2], default=2, help=axis_help)
enc_parser.add_argument('--sort-rows', help='Sort rows', action='store_true')
enc_parser.add_argument('--sort-cols', help='Sort columns', action='store_true')
# enc_parser.add_argument('-t', '--transpose', help='Tranpose binary matrix', action='store_true')
enc_parser.add_argument('--dist', help='Distance to compute cost matrix', choices=AVAIL_DIST, default='ham')
enc_parser.add_argument('--solver', help='Traveling Salesman Problem solver', choices=AVAIL_SOLVERS, default='nn')
enc_parser.add_argument('--preset-mode', type=int, choices=[0, 1, 2], default=0, help=preset_mode_help)
enc_parser.add_argument('--encoder', choices=AVAIL_CODECS, default='jbig', help=encoder_help)
enc_parser.add_argument('input', metavar='input_file_path', help='Input file path.')
enc_parser.add_argument('output', metavar='output_file_path', help='Output file path.')

#? Decode
dec_parser = subparsers.add_parser('decode')
dec_parser.add_argument('--pos', type=int, help='The start and end position of genomic variants to be decoded.', nargs=2)
dec_parser.add_argument('--samples', type=str, help='A semicolon separated list of samples to be decoded. Default: all samples.')
dec_parser.add_argument('input', metavar='input_file_path', help='Input file path.')
dec_parser.add_argument('output', metavar='output_file_path', help='Output file path.', nargs='?')

# #? Random Access
# ra_parser = subparsers.add_parser('random-access')
# ra_parser.add_argument('input', metavar='input_file_path', help='input file path')
# ra_parser.add_argument('--pos', type=int, help='The start and end position of genomic variants to be decoded.', nargs=2)
# ra_parser.add_argument('--samples', type=str, help='List of samples to be decoded. Default: all samples', nargs='+')

args = parser.parse_args()

gvc.app.run(args, run_as_module=False)
