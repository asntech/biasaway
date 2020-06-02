#!/usr/bin/env python

"""
This is a setup script for BiasAway: a tool for DNA sequence background generation

This code is free software; you can redistribute it and/or modify it under the terms of the
BSD License (see the file LICENSE.md included with the distribution).

@author: Aziz Khan
@email: aziz.khan@ncmm.uio.no
"""
import os
from distutils.core import setup
from setuptools import find_packages
from biasaway import __version__ as VERSION


CLASSIFIERS = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

install_requires = [
    'biopython',
    'numpy',
]

#def readme():
#    with open('README.rst') as f:
#        return f.read()

def readme(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()

setup(
    name="biasaway",
    description="a tool to generate nucleotide composition-matched DNA sequences",
    version=VERSION,
    author="Anthony Mathelier and Aziz Khan",
    license='GPL',
    platforms='linux/unix',
    author_email="anthony.mathelier@ncmm.uio.no",
    url="https://bitbucket.org/CBGR/biasaway/src/master/",
    long_description=readme("README.rst"),
    long_description_content_type='text/markdown',
    package_dir={'biasaway': 'biasaway'},

    packages=['biasaway',
        #'biasaway.example_data'
        ],

    scripts=['biasaway/biasaway',
                   ],
    include_package_data=True,
    install_requires = install_requires,
    classifiers=CLASSIFIERS,
)
