#!/usr/bin/env python

# --------------------------------------------------------------- #
# Phylogenetic Clustering by Linear Integer Programming (PhyCLIP) #
# Authors: Alvin X. Han and Edyth Parker                          #
# --------------------------------------------------------------- #

import os
from setuptools import setup

setup(
    name = 'PhyCLIP',
    version = '1.1',
    author = 'Alvin X. Han and Edyth Parker',
    author_email = 'hanxc@bii.a-star.edu.sg',
    description = ('Phylogenetic clustering by linear integer programming.'),
    keywords = 'Phylogenetics, Clustering',
    url = 'https://github.com/alvinxhan/PhyCLIP',
    packages = ['phyclip_modules'],
    scripts = ['bin/phyclip.py']
)
