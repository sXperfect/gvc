from enum import IntEnum

class DataUnitType(IntEnum):
    PARAMETER_SET = 0
    ACCESS_UNIT = 1

class ContentID(IntEnum):
    GENOTYPE = 0

class ConcatAxis(IntEnum):
    ROW = 0
    COL = 1
    NONE = 2
    TOTAL = 3

class BinarizationFlag(IntEnum):
    pass
