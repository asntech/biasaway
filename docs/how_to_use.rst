====================
How to use BiasAway
====================

Once you have installed BiasAway, you can type:

.. code-block:: bash

    biasaway --help

It will print the main help, which lists the six subcommands/modules: ``m``, ``f``, ``d``, ``w``, ``g``, and ``c``.

.. code-block:: bash

	usage: biasaway <subcommand> [options]

	positional arguments <subcommand>: {m,f,d,w,g,c}

		List of subcommands
		m 	mononucleotide shuffling
		d 	dinucleotide shuffling
		f 	mononucleotide shuffling within a sliding window
		w 	dinucleotide shuffling within a sliding window
		g 	mononucleotide distribution matched
		c 	mononucleotide distribution within a sliding window matched

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit

To view the help for the individual subcommands, please type:

.. note:: Please check ``BiasAway modules`` to see a detailed summary of available **options**.

To view ``m`` module help, type

.. code-block:: bash

	biasaway m --help

To view ``f`` module help, type

.. code-block:: bash

	biasaway f --help

To view ``d`` module help, type

.. code-block:: bash

	biasaway d --help

To view ``w`` module help, type

.. code-block:: bash

	biasaway w --help

To view ``g`` module help, type

.. code-block:: bash

	biasaway g --help


To view ``c`` module help, type

.. code-block:: bash

	biasaway c --help
