import inspect
import os
import logging as log

_self_filename = inspect.getfile(inspect.currentframe())

PROGRAM_NAME = 'gvc'
PROGRAM_DESC = 'Genomic Variant Codec'
THIRD_PARTY_DIR = os.path.join(os.path.dirname(os.path.dirname(_self_filename)), 'third_party')

LOG_LEVELS = {
    'critical': log.CRITICAL,
    'error': log.ERROR,
    'warning': log.WARNING,
    'info': log.INFO,
    'debug': log.DEBUG,
}
