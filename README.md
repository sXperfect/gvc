# Genomic Variant Codec (GVC)

[![CI](https://github.com/sXperfect/gvc/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/sXperfect/gvc/actions/workflows/ci.yml)

Open Source Genotype Compressor

## Usage policy
---

The open source GVC codec is made available before scientific publication.

This pre-publication software is preliminary and may contain errors.
The software is provided in good faith, but without any express or implied warranties.
We refer the reader to our [license](LICENSE).

The goal of our policy is that early release should enable the progress of science.
We kindly ask to refrain from publishing analyses that were conducted using this software while its development is in progress.

## Dependencies
---

GVC 1.0.x supports Python 3.8 or newer. CMake and a C/C++ compiler are
required for the native components. Future minor release lines may raise the
minimum Python version; older release lines remain available for legacy Python
environments. See [docs/design/versioning.md](docs/design/versioning.md).
For anaconda or conda user, CMAKE, gcc and gxx libraries are required and can be installed through: `conda install -c conda-forge cmake gxx_linux-64 gcc_linux-64`.
The core numerical dependencies are declared in `pyproject.toml`. Install
optional integrations explicitly:

```bash
python -m pip install -e ".[vcf]"      # VCF/BCF parsing
python -m pip install -e ".[speed]"    # Numba acceleration
python -m pip install -e ".[all]"      # all optional runtime integrations
python -m pip install -e ".[test]"     # test/build tooling
```

`requirements.txt` remains a full-feature compatibility install for legacy
workflows.

## Building
---

Clone this repository:

    git clone https://github.com/sXperfect/gvc

For a development installation:

    python -m pip install -e .

Build all native components and run the local verification gate with:

    bash setup.sh
    ./scripts/verify.sh

The local verification gate compiles the Cython extensions and the standalone
CMake helper before running the test suite.

### Entropy Codec
---

In order to encode or decode the payloads based on JBIG codec, an external executable is required.
You can use any of the existing and publicly available JBIG codec implementation.
We provide an example on how to integrate JBIG-based codec [here](JBIG.md).

Generic compressors, such as LZMA or BZIP2, are supported.
Please refer to this [documentation](CODEC.md) for integration.

## Usage

<!-- In order to use GVC, you should activate the virtual environment first.
To activate the environment, we can use the following command in the root folder: `source tmp/env/bin/activate`.
Alternatively, you can replace `python3` in the following commands with `tmp/venv/bin/python3`. -->
<!-- If you install the dependencies manually, you do not have to activate the environment. -->

---
Compress a VCF file with default options (an example VCF file can be found in the `tests` folder): 
```
python3 -m gvc encode variant_calls.vcf compressed_genotypes.gvc
```

A list of options can be obtained via:
```
python3 -m gvc encode --help
```

Decode a compressed VCF file with default options: 
```
python3 -m gvc decode compressed_genotypes.gvc decoded_genotypes.txt
```

A list of options can be obtained via:
```
python3 -m gvc decode --help
```

For random access to a subset of compressed genotypes, additional options must be passed to the `python3 -m gvc decode` command:

```
python3 -m gvc decode --pos 1 10 --sample SAMPLE01 compressed_genotypes.gvc decoded_genotypes.txt
```
