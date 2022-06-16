# Genomic Variant Codec (GVC)

Open Source Genotype Compressor

---

## Usage policy

The open source GVC codec is made available before scientific publication.

This pre-publication software is preliminary and may contain errors.
The software is provided in good faith, but without any express or implied warranties.
We refer the reader to our [license](LICENSE).

The goal of our policy is that early release should enable the progress of science.
We kindly ask to refrain from publishing analyses that were conducted using this software while its development is in progress.

--- 
## Dependencies

* Python>=3.6
* virtualenv
* numpy
* cython
* numba
* cyvcf2
* scipy
* tspsolve
* Pillow

---
## Building

Clone this repository:

    git clone https://github.com/sXperfect/gvc

Run setup script `setup.sh`

    bash setup.sh

This will create automaticaly a virtual environment with all dependencies installed.
The newly created virtual environment will be located in `tmp/venv`.

---
### Entropy Codec

In order to encode or decode the payloads based on JBIG codec, an external executable is required.
You can use any of the existing and publicly available JBIG codec implementation.
We provide an example on how to integrate JBIG-based codec [here](JBIG).

Assume that a generic compressor should be used, such as LZMA or BZIP2, please refer to the documentation [here](CODEC)