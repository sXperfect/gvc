
In order to encode or decode the payloads based on JBIG codec, an external executable is required.
You can use any of the existing and publicly available JBIG codec implementation.
Here we provides an example based on [JBIG-KIT](https://www.cl.cam.ac.uk/~mgk25/jbigkit/).
First, jbig must be downloaded into folder `third_party` and compiled.

In order to integrate jbig-kit, two functions are required.
One for the encoding process and the other one for the decoding process.
The functions call jbig executeable and then the compressed payload is read.
We create a new file called `jbigkit.py`

```python
import os
import logging as log
import subprocess as sp
import tempfile as tmp

from PIL import Image
import numpy as np

from .jbig import get_shape
from .. import settings
from .. import utils

# Allow super large image decompression. DO NOT REMOVE!
Image.MAX_IMAGE_PIXELS = np.inf

JBIGKIT_PATH = os.path.join(settings.THIRD_PARTY_DIR, 'jbigkit-2.1', 'pbmtools')
ENCODER_FPATH = os.path.join(JBIGKIT_PATH, "pbmtojbg85") #? Specify path to encoder
DECODER_FPATH = os.path.join(JBIGKIT_PATH, "jbgtopbm85") #? Specify path to decoder
utils.check_executable(ENCODER_FPATH)
utils.check_executable(DECODER_FPATH)

JBG_LRLTWO = 0x40
JBG_VLENGTH = 0x20
JBG_TPBON = 0x08

def encode(
    matrix:np.ndarray
) -> bytes:
    """
    Encode a binary matrix as JBIG1 bytestream.

    Parameters
    ----------
    matrix : ndarray (2d)
        binary_matrix: The binary matrix to encode. Must contain integers equal to 0 or 1

    Returns
    -------
    jbig1_binary: binary
        The JBIG1 bytestream.

    """
    with tmp.TemporaryDirectory(prefix='PBM2JBG_') as dpath:
        
        pbm_fpath = os.path.join(dpath, 'tmp.pbm')
        jbg_fpath = os.path.join(dpath, 'tmp.jbg')

        # Store matrix as PBM file
        pbm_im = Image.fromarray(matrix.astype(bool))
        pbm_im.save(pbm_fpath)

        flags = JBG_TPBON
        l0 = (1<<32)-1

        args = [
            "-p", str(flags),
            "-s", str(l0)
        ]

        proc = sp.run([ENCODER_FPATH, *args, pbm_fpath, jbg_fpath], stdout=sp.PIPE, stderr=sp.PIPE)

        if proc.returncode != 0:
            log.error(proc.stderr)
            raise RuntimeError(proc.stderr.decode("UTF-8"))

        with open(jbg_fpath, 'rb') as jbig_f:
            jbig_payload = jbig_f.read()

    return jbig_payload

def decode(
    jbig_payload:bytes
) -> np.ndarray:

    """
    Decode a JBIG1 bytestream as binary matrix.

    Parameters
    ----------
    jbig_payloadtream: bytes
        The JBIG1 bytestream.

    Returns
    -------
    bin_matrix: np.ndarray
        The binary matrix
    """

    with tmp.TemporaryDirectory(prefix='JBG2PBM_') as dpath:
        jbg_fpath = os.path.join(dpath, 'tmp.jbg')
        pbm_fpath = os.path.join(dpath, 'tmp.pbm')

        with open(jbg_fpath, 'wb') as jbig1_f:
            jbig1_f.write(jbig_payload)

        __, ncols = get_shape(jbig_payload[:12])

        args = [
            "-x", str(ncols),
            "-B", str(2**30),
        ]

        proc = sp.run([DECODER_FPATH, *args, jbg_fpath, pbm_fpath], stdout=sp.PIPE, stderr=sp.PIPE)

        if proc.returncode != 0:
            log.error(proc.stderr)
            raise RuntimeError(proc.stderr.decode("UTF-8"))

        pbm_im = Image.open(pbm_fpath)

    return np.asarray(pbm_im)
```
in `gvc/codec` folder. Both functions should specify all necessary flags. 
Internally, GVC does not extract the header of the payload, thus JBIG parameters are independent from GVC.

We then register the encoding and decoding functions in `__init__.py` file located in `gvc/codec`:

```python

from . import jbigkit #? Import jbigkit

#? If a new codec is added, please update data_structure.consts too
MAT_CODECS = {
    CodecID.JBIG1 : {
        "name": "jbig",
        "encoder": jbigkit.encode, #? Add encode function
        "decoder": jbigkit.decode, #? Add decode function
    }
}
```

Now JBIG is ready to use.