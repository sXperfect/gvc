# Genomic Variant Codec (GVC)

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

Both python version 3.7 or newer and python virtual environment are required.
See [requirements.txt](requirements.txt) for the list of required python libraries.
For python version 3.6 or lower, an additional python package `dataclass` is required.

## Building
---

Clone this repository:

    git clone https://github.com/sXperfect/gvc

Run setup script `setup.sh`

    bash setup.sh

This will create automaticaly a virtual environment with all dependencies installed located in `tmp/venv`.

### Entropy Codec
---

In order to encode or decode the payloads based on JBIG codec, an external executable is required.
You can use any of the existing and publicly available JBIG codec implementation.
We provide an example on how to integrate JBIG-based codec [here](JBIG.md).

Generic compressors, such as LZMA or BZIP2, are supported.
Please refer to this [documentation](CODEC.md) for integration.

## Usage

In order to use GVC, you should activate the virtual environment first.
To activate the environment, we can use the following command in the root folder: `source tmp/env/bin/activate`.
Alternatively, you can replace `python3` in the following commands with `tmp/venv/bin/python3`.

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
