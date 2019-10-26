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
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="biasaway",
    description="a tool to generate composition-matched background sequence sets",
    version=VERSION,
    author="Aziz Khan",
    license='MIT',
    platforms='linux/unix',
    author_email="azez.khan@gmail.com",
    url="https://github.com/asntech/biasaway",
    long_description=readme("README.rst"),
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
