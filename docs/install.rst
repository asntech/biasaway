============
Installation
============
BiasAway is available on `PyPi <https://pypi.python.org/pypi/biasaway>`_,
through `Bioconda <https://bioconda.github.io/recipes/biasaway/README.html>`_,
and the source code is available on `bitbucket
<https://bitbucket.org/cbgr/biasaway>`_. BiasAway takes care of the
installation of all the required python modules. If you already have a working
installation of python, the easiest way to install the required python modules
is by installing biasaway using ``pip``. 

If you are setting up Python for the first time, we recommend to install it
using the `Conda or Miniconda Python distribution
<https://conda.io/docs/user-guide/install/index.html>`_. This comes with
several helpful scientific and data processing libraries available for
platforms including Windows, Mac OSX, and Linux.

You can use one of the following ways to install BiasAway.

Quick installation
==================

Prerequisites
=============
BiasAway requires the following Python modules:

	* biopython: https://biopython.org
	* numpy: https://numpy.org
    * matplotlib: https://matplotlib.org/
    * seaborn: https://seaborn.pydata.org/

Install biopython, numpy, matplotlib, and seaborn
-------------------------------------------------
BiasAway uses `biopython <https://biopython.org>`_, `numpy
<https://numpy.org>`_, `matplotlib <https://matplotlib.org>`_, and `seaborn
<https://seaborn.pydata.org>`_ you can install them using `pip` or `conda`.

.. code-block:: bash
    conda install biopython numpy matplotlib seaborn

.. note:: If you install using ``pip`` or ``bioconda`` prerequisites will be installed. 


Install BiasAway using `conda`
==============================
BiasAway is available on `Bioconda <https://anaconda.org/bioconda/biasaway>`_ for installation via ``conda``.

.. code-block:: bash

	conda install -c bioconda biasaway


Install BiasAway using `pip`
============================
BiasAway is available on `PyPi <https://pypi.org/project/biasaway/>`_ for installation via ``pip``.

.. code-block:: bash

	pip install biasaway


Install BiasAway from source
=============================
You can install the development version by using ``git`` from our bitbucket
repository at https://bitbucket.org/CBGR/biasaway. 


Install development version from `Bitbucket`
--------------------------------------------

If you have `git` installed, use this:

.. code-block:: bash

    git clone https://bitbucket.org/CBGR/biasaway.git
    cd biasaway
    python setup.py sdist install
