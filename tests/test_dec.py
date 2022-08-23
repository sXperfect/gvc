            "0/0\t0|1\t3/0",
            "3|1\t0/2\t1|1",
            "3|1\t0/2\t1|2",
            "0/0\t0|1\t3/1",
            "0/0\t.|.\t3/1",
            
from os import getcwd
from os.path import join
import itertools as it
import numpy as np
from gvc.common import create_parameter_set
from gvc.sort import sort
from gvc import reader
from gvc import data_structures as ds
from gvc.codec import encode, decode
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.encoder import run_core, binarize_allele_matrix
from gvc.decoder import decode_encoded_variants
import unittest


class TestDecode(unittest.TestCase):