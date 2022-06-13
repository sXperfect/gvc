from setuptools import setup
# from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

setup(
    name='cdebinarize',
    ext_modules=cythonize(
        "gvc/cdebinarize.pyx",
        language_level="3",
        # language="c++",
        annotate=True,
    ),
    # ext_modules = cythonize(extensions),
    zip_safe=False,
    include_dirs=[numpy.get_include()]
)

setup(
    name='c_gpbwt',
    ext_modules=cythonize(
        "gvc/c_pbwt.pyx",
        language_level="3",
        # language="c++",
        annotate=True,
    ),
    # ext_modules = cythonize(extensions),
    zip_safe=False,
    include_dirs=[numpy.get_include()]
)

setup(
    name='c_rle',
    ext_modules=cythonize(
        "gvc/c_rle.pyx",
        language_level="3",
        # language="c++",
        annotate=True,
    ),
    # ext_modules = cythonize(extensions),
    zip_safe=False,
    include_dirs=[numpy.get_include()]
)