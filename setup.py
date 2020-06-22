#!/usr/bin/env python

"""
This is a setup script for BiasAway: a tool for DNA sequence background
generation

This code is free software; you can redistribute it and/or modify it under the
terms of the BSD License (see the file LICENSE included with the distribution).

@author: Anthony Mathelier
@email: anthony.mathelier@ncmm.uio.no
"""
import os
from setuptools import setup
from biasaway import __version__ as VERSION


CLASSIFIERS = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

install_requires = [
    'biopython',
    'numpy',
    'matplotlib',
    'seaborn',
    'ushuffle',
    'sklearn',
    'scipy'
]


def readme(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()


descr = "a tool to generate nucleotide composition-matched DNA sequences"
setup(
    name="biasaway",
    description=descr,
    version=VERSION,
    author="Anthony Mathelier and Aziz Khan",
    license='GPL',
    platforms='linux/unix',
    author_email="anthony.mathelier@ncmm.uio.no",
    url="https://bitbucket.org/CBGR/biasaway/src/master/",
    long_description=readme("README.rst"),
    long_description_content_type='text/x-rst',
    package_dir={'biasaway': 'biasaway'},
    packages=['biasaway'],
    scripts=['biasaway/biasaway'],
    include_package_data=True,
    install_requires=install_requires,
    classifiers=CLASSIFIERS,
)
