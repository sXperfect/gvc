from os.path import join, basename
import glob
import shutil
from setuptools import setup
import numpy as np
from distutils.extension import Extension
from Cython.Build import build_ext

fname = "cdebinarize"
setup(
    name=fname,
    ext_modules=[
        Extension(
            fname,
            sources=["gvc/" + fname + ".pyx"],
            extra_compile_args=['-fopenmp', '-O3', '-march=native'],
            extra_link_args=['-fopenmp', '-O3', '-march=native'],
            language='c++'
        )
    ],
    cmdclass = {'build_ext': build_ext},
    include_dirs=[np.get_include()],
    zip_safe=False,
)


fname = "cquery"
setup(
    name=fname,
    ext_modules=[
        Extension(
            fname,
            sources=["gvc/" + fname + ".pyx"],
            extra_compile_args=['-fopenmp', '-O3', '-march=native'],
            extra_link_args=['-fopenmp', '-O3', '-march=native'],
            language='c++'
        )
    ],
    cmdclass = {'build_ext': build_ext},
    include_dirs=[np.get_include()],
    zip_safe=False,
)

fpath = glob.glob(fname + '*.so')[0]
so_fname = basename(fpath)
shutil.move(fpath, join('gvc', so_fname))
