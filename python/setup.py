import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "stokesdlpy"

list_files=[]
list_files.append('../stokes_sherman_lauricella_solver.f')

FLIBS = os.getenv('STLIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None,FLIBS))
FLIBS.append('../lib-static/libstokes2ddl.a')
if platform == "darwin":
    FLIBS.append('-L/usr/local/lib')

stok = []
stok.append('get_bdry_data')
stok.append('stokes_gmres')
stok.append('eval_vel')


ext_stok = Extension(
    name='stokesdl',
    sources=list_files,
    f2py_options=['only:']+stok+[':'],
#    extra_f77_compile_args=FFLAGS,
    extra_f90_compile_args=["-std=legacy"],
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Manas Rachh",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for stokes 2d dl",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_stok],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
