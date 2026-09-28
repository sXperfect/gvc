from __future__ import annotations

import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

compile_args = ["/O2"] if os.name == "nt" else ["-O3"]

extensions = [
    Extension(
        "gvc.data_structures.crc_id",
        ["gvc/data_structures/crc_id.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=compile_args,
        language="c++",
    ),
    Extension(
        "gvc.cdebinarize",
        ["gvc/cdebinarize.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=compile_args,
        language="c++",
    ),
    Extension(
        "gvc.cquery",
        ["gvc/cquery.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=compile_args,
        language="c++",
    ),
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),
    zip_safe=False,
)
