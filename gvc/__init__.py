"""Genomic Variant Codec."""

from importlib import metadata

try:
    __version__ = metadata.version("gvc")
except metadata.PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = ["__version__"]
