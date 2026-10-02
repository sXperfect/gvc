# JBIG-KIT integration

GVC's historical `CodecID.JBIG1` payload uses an external JBIG1/T.85
implementation. GVC 1.0.x now contains a maintained subprocess integration in
`gvc.codec.jbigkit`; users no longer need to copy the old example module into the
package.

## Install JBIG-KIT

The integration expects these executables:

- `pbmtojbg85` for encoding;
- `jbgtopbm85` for decoding.

On Ubuntu/Debian systems that provide JBIG-KIT binaries, install the
`jbigkit-bin` package. Alternatively, build JBIG-KIT 2.1 and point GVC at the
resulting executables.

Pillow is used for PBM interchange and is available through the GVC JBIG/all
optional dependency groups:

```bash
python -m pip install ".[jbig]"
```

## Executable discovery

By default GVC searches `PATH` for `pbmtojbg85` and `jbgtopbm85`.

Custom locations can be configured without editing GVC:

```bash
export GVC_JBIG_ENCODER=/opt/jbigkit/pbmtools/pbmtojbg85
export GVC_JBIG_DECODER=/opt/jbigkit/pbmtools/jbgtopbm85
```

The subprocess timeout defaults to 60 seconds per matrix and can be changed
with:

```bash
export GVC_JBIG_TIMEOUT=300
```

A non-positive timeout is rejected.

## Python API

The codec is registered automatically as GVC's `jbig` codec. Direct use is also
available:

```python
import numpy as np
from gvc.codec import jbigkit

matrix = np.array([[0, 1], [1, 0]], dtype=bool)
payload = jbigkit.encode(matrix)
restored = jbigkit.decode(payload)
```

The wrapper validates:

- two-dimensional binary input;
- executable existence and execute permission;
- subprocess timeout and non-zero exits;
- non-truncated JBIG headers;
- encoded/decoded matrix shape consistency;
- readable decoder PBM output.

External command failures raise `JBIGKitError` with the command status and stderr.

## Multiprocessing

The maintained functions are module-level and therefore compatible with
GVC's spawn/fork multiprocessing architecture. They do not require inherited
runtime registry mutations.

If an application replaces the codec dynamically, spawn/forkserver workers can
recreate that state with `multiprocessing_initializer`; see
`docs/design/multiprocessing.md`.

## Historical compatibility verification

The original LUH repository used `tests/test_block01.vcf.gz` as its main
large encode/decode fixture. It is intentionally not committed again to this
repository because it is approximately 7.5 MB compressed.

The Python 3.8 release CI fetches the fixture from the pinned original commit:

```text
tnt-LUH/gvc
commit: f9af2127a2ff0b87727924860e33f1905fd507cd
blob:   af45a419e46563906ac51fad0be869291316cea9
```

For a larger offline verification, download that fixture and run:

```bash
python scripts/verify_historical.py tests/test_block01.vcf.gz --max-blocks 0
```

To additionally reproduce the original row/column sorting combinations:

```bash
python scripts/verify_historical.py tests/test_block01.vcf.gz \
  --max-blocks 0 \
  --include-sorting
```

These full runs can be substantially more expensive than the bounded hosted
release gate.

## Format behavior

GVC stores the JBIG payload bytes rather than normalizing external codec
parameters. The serialized GVC framing therefore remains compatible with the
historical design while the executable integration, validation, and failure
handling are maintained by GVC.
