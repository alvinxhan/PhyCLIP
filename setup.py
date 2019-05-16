#!/usr/bin/env python

# --------------------------------------------------------------- #
# Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) #
# Authors: Alvin X. Han and Edyth Parker                          #
# --------------------------------------------------------------- #

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


ext_modules=[
    Extension("phyilpx",
              ["phyilpx/*.pyx"],
              include_dirs=[np.get_include()],
              extra_link_args=['-L/usr/lib/x86_64-linux-gnu/']
              )
]

setup(
    name = 'PhyCLIP',
    version = '2.0',
    author = 'Alvin X. Han and Edyth Parker',
    author_email = 'hanxc@bii.a-star.edu.sg',
    description = ('Phylogenetic clustering by linear integer programming.'),
    keywords = 'Phylogenetics, Clustering',
    url = 'https://github.com/alvinxhan/PhyCLIP',
    include_dirs=[np.get_include()],
    ext_modules=cythonize("phyilpx/*.pyx"),
    packages = ['phyclip_modulex'],
    scripts = ['bin/phyclip.py']
)
