# Generic Compressor

Here we provide an example based on LZMA.
First, we add a new constant calld `LZMA` to the `CodecID` in `gvc/data_structures/consts.py` class to accomodate the new codec.

```python
class CodecID(IntEnum):
    JBIG1 = 0
    GABAC = 1
    LZMA = 2
```

Similar to [JBIG](JBIG) implementation, both encode and decode functions must be provided.
Furthermore, matrix serialization and additional information is required such as 
the number of rows and the number of rows in order to reconstruct the serialized
matrix.

We then register the encode and decode functions in `__init__.py` file located in `gvc/codec`.
This is done by adding a new entries in `MAT_CODECS` variable:

```python
from . import lzma #? Import LZMA

#? If a new codec is added, please update data_structure.consts too
MAT_CODECS = {
    CodecID.LZMA : { #? Add new dict entries for additional codec
        "name": "LZMA",
        "encoder": lzma.encode, #? Add encode function
        "decoder": lzma.decode, #? Add decode function
    }
}
```

