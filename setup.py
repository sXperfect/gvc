import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup


def _extension(name, source):
    compile_args = ["-O3"]
    link_args = []
    if sys.platform.startswith("linux"):
        compile_args.append("-fopenmp")
        link_args.append("-fopenmp")

    return Extension(
        name,
        sources=[source],
        include_dirs=[np.get_include()],
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        language="c++",
    )


extensions = [
    _extension("gvc.data_structures.crc_id", "gvc/data_structures/crc_id.pyx"),
    _extension("gvc.cdebinarize", "gvc/cdebinarize.pyx"),
    _extension("gvc.cquery", "gvc/cquery.pyx"),
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),
)
