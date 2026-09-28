# Genomic Variant Codec (GVC)

[![CI](https://github.com/sXperfect/gvc/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/sXperfect/gvc/actions/workflows/ci.yml)

GVC is a research codec for compact representation and random-access processing of genotype data. This repository modernizes the historical GVC implementation while preserving its BSD license and codec/data-structure lineage.

## Release and Python support

The current modernization line is **GVC 1.x**, with a minimum supported Python version of **3.8**.

| GVC line | Minimum Python | Status |
| --- | ---: | --- |
| 1.x | 3.8 | Active compatibility/modernization line |
| 2.x | Higher floor, to be decided | Future major modernization |

Raising the minimum Python version is treated as a **major-release change**. A future GVC 2.x may adopt a newer Python baseline without forcing users of Python 3.8 to upgrade away from the final compatible 1.x release.

See [Python support policy](docs/development/python-support.md).

## Modernization status

The 1.x work focuses on making the historical implementation reproducible and testable before changing its compatibility floor:

- PEP 517/518 packaging with explicit Python metadata;
- self-contained pytest coverage for bitstreams, binarization, sorting, and codec round trips;
- removal of obsolete mandatory dependencies where a small maintained implementation is sufficient;
- optional VCF and acceleration dependencies;
- portable Cython/CMake build paths;
- a local-first CI gate modeled after the gunz-utils CI design;
- explicit preservation of historical licensing and provenance.

The on-disk format and core algorithms remain compatibility-sensitive. Algorithmic changes should be covered by regression tests before they are adopted.

## Installation

Core development installation:

    python -m pip install -e ".[test]"

VCF parsing support is optional:

    python -m pip install -e ".[test,vcf]"

Optional Numba acceleration:

    python -m pip install -e ".[test,speed]"

A C/C++ compiler is required to build the Cython extensions. CMake is required only for the standalone helper under library/libgvc.

## Verification

Run the same consolidated gate used by hosted CI:

    bash scripts/verify.sh

Focused checks are available through:

    python scripts/ci.py syntax
    python scripts/ci.py native
    python scripts/ci.py test
    python scripts/ci.py all

## Command line

Show the command-line interface:

    python -m gvc --help

Encoding and decoding retain the historical interface:

    python -m gvc encode input.vcf output.gvc
    python -m gvc decode input.gvc output.txt

VCF-backed commands require the vcf extra.

## Entropy-codec integration

The historical GVC codec registry defines the JBIG codec contract but does not bundle a JBIG executable. Production JBIG encode/decode therefore still requires an integration that supplies the codec functions described in [JBIG.md](JBIG.md).

The unit suite deliberately does **not** depend on an external JBIG binary. It injects a deterministic lossless test codec so that GVC's own binarization, sorting, payload, and reconstruction logic can be tested independently of third-party codec installation.

## Repository provenance

This repository is a modernization of the historical open-source GVC implementation previously maintained at tnt-LUH/gvc. The original BSD license and copyright notices are retained in [LICENSE](LICENSE), and the historical README is preserved at [docs/history/README-upstream.md](docs/history/README-upstream.md).

See [NOTICE.md](NOTICE.md) for provenance details.
