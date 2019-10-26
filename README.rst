BiasAway
--------

	a tool to generate composition-matched background sequence sets

.. image:: https://travis-ci.org/asntech/biasaway.svg?branch=master
    :target: https://travis-ci.org/asntech/biasaway

.. image:: https://img.shields.io/pypi/pyversions/biasaway.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/biasaway.svg
    :target: https://pypi.python.org/pypi/biasaway

.. image:: https://anaconda.org/bioconda/biasaway/badges/version.svg
	:target: https://anaconda.org/bioconda/biasaway

.. image:: https://anaconda.org/bioconda/biasaway/badges/downloads.svg
    :target: https://bioconda.github.io/recipes/biasaway/README.html

.. image:: https://anaconda.org/bioconda/biasaway/badges/installer/conda.svg
	:target: https://conda.anaconda.org/bioconda

.. image:: https://img.shields.io/github/issues/asntech/biasaway.svg
	:target: https://github.com/asntech/biasaway/issues


Documentation
=============

**A detailed documentation is available in different formats:**  `HTML <http://biasaway.readthedocs.org>`_ | `PDF <http://readthedocs.org/projects/biasaway/downloads/pdf/latest/>`_ | `ePUB <http://readthedocs.org/projects/biasaway/downloads/epub/latest/>`_


Installation
============

Quick installation using Conda
------------------------------

.. code-block:: bash

	conda install -c bioconda biasaway

This will install all the dependencies and you are ready to use BiasAway.

Install using `pip`
-------------------
You can install BiasAway from PyPi using pip.

Install from PyPi::

	pip install biasaway

Note: If you install using pip, make sure to install BEDTools and R packages listed below. 

BiasAway requires the following Python modules and R packages:

	* Python (v2.7): https://www.python.org
	* Biopython: https://biopython.org
	* Numpy (<=v1.16.5): https://numpy.org

A quick installation using ``pip``, you can also use conda.

.. code-block:: bash

    pip install biopython

.. code-block:: bash

    pip install numpy


Install BiasAway from source
=============================
You can install a development version by using ``git`` from GitHub or Bitbucket.


Install development version from `Bitbucket`
--------------------------------------------

If you have `git` installed, use this:

.. code-block:: bash

    git clone https://bitbucket.org/CBGR/biasaway.git
    cd biasaway
    python setup.py sdist install

Install development version from `GitHub`
-----------------------------------------
If you have `git` installed, use this:

.. code-block:: bash

    git clone https://github.com/asntech/biasaway.git
    cd biasaway
    python setup.py sdist install

How to use BiasAway
====================
Once you have installed biasaway, you can type:

.. code-block:: bash

	biasaway --help

This will show the main help, which lists the three subcommands/modules: ``m``, ``f``, ``d``, ``w``, ``g``, and ``c``.

.. code-block:: bash

	usage: biasaway <subcommand> [options]

		positional arguments <subcommand>: {m,f,d,w,g,c}

		List of subcommands
		m 	mono-nucleotide shuffling generator
		f 	mono-nucleotide shuffling within a sliding window generator
		d 	di-nucleotide shuffling generator
		w 	di-nucleotide shuffling within a sliding window generator
		g 	%GC distribution-based background chooser
		c 	GC distribution and %GC composition within a sliding window background chooser

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit


to see the help for the three subcommands ``m``, ``f``, ``d``, ``w``, ``g``, and ``c`` type:

.. code-block:: bash
	
	biasaway m --help

	biasaway f --help

	biasaway d --help

	biasaway w --help

	biasaway g --help

	biasaway c --help


Interactive Web App
=====================
BiasAway Web App is freely available at: 

The source code for the web app is available at https://github.com/asntech/BiasAwayApp

Support
========
If you have questions, or found any bug in the program, please write to us at ``aziz.khan[at]ncmm.uio.no``

