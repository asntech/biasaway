============
Installation
============
BiasAway is available on `PyPi <https://pypi.python.org/pypi/biasaway>`_, through `Bioconda <https://bioconda.github.io/recipes/biasaway/README.html>`_, and source code available on `GitHub <https://github.com/asntech/biasaway>`_ and `Bitbucket <https://bitbucket.org/CBGR/biasaway>`_. BiasAway takes care of the installation of all the required Python modules. If you already have a working installation of Python, the easiest way to install the required Python modules is by installing BiasAway using ``pip``. 

If you're setting up Python for the first time, we recommend to install it using the `Conda or Miniconda Python distribution <https://conda.io/docs/user-guide/install/index.html>`_. This comes with several helpful scientific and data processing libraries, and available for platforms including Windows, Mac OSX and Linux.

You can use one of the following ways to install BiasAway.

Quick installation
==================

Prerequisites
=============
BiasAway requires the following Python modules:

	* Biopython: https://biopython.org
	* Numpy: https://numpy.org

Install Biopython and Numpy
----------------------------
BiasAway is using `biopython <https://biopython.org>`_ and `numpy <https://numpy.org>`_, you can install using `pip`.

.. note:: If you install using ``pip`` or ``bioconda`` prerequisites will be installed. 

.. code-block:: bash

	pip install biopython

.. code-block:: bash

	pip install numpy

Install uisng Conda
--------------------
We highly recommend to install BiasAway using Conda, this will take care of the dependencies. If you already have Conda or Miniconda installed, go ahead and use the below command.

.. code-block:: bash

	conda install -c bioconda biasaway

.. note:: This will install all the dependencies and you are ready to use **BiasAway**.

Install using `pip`
-------------------
You can install BiasAway from PyPi using pip.

.. code-block:: bash

	pip install biasaway


Install BiasAway from source
=============================
You can install a development version by using ``git`` from our bitbucket repository at https://bitbucket.org/CBGR/biasaway or Github. 


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