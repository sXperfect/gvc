from os import remove
from os.path import join, basename, dirname
import glob
import shutil
from setuptools import setup
import numpy as np
from distutils.extension import Extension
from Cython.Build import build_ext

fpaths = [
    "gvc/data_structures/crc_id.pyx",
    "gvc/cdebinarize.pyx",
    "gvc/cquery.pyx"
]


for fpath in fpaths:
    fname = basename(fpath).replace('.pyx', '')
    
    setup(
        name=fname,
        ext_modules=[
            Extension(
                fname,
                sources=[fpath],
                extra_compile_args=['-fopenmp', '-O3', '-march=native'],
                extra_link_args=['-fopenmp', '-O3', '-march=native'],
                language='c++'
            )
        ],
        cmdclass = {'build_ext': build_ext},
        include_dirs=[np.get_include()],
        zip_safe=False,
    )

    so_fpath = glob.glob(fname + '*.so')[0]
    target_dir = dirname(fpath)
    
    try:
        remove(join(target_dir, basename(so_fpath)))
    except:
        pass
    shutil.move(so_fpath, target_dir)
