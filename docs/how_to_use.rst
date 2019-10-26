====================
How to use BiasAway
====================

Once you have installed BiasAway, you can type:

.. code-block:: bash

    biasaway --help

This will show the main help, which lists the three subcommands/modules: ``m``, ``f``, ``d``, ``w``, ``g``, and ``c``.

.. code-block:: bash

	usage: biasaway <subcommand> [options]
	    
	positional arguments <subcommand>:
		{m,f,d,w,g,c}
			List of subcommands
	m		mono-nucleotide shuffling generator
	f		mono-nucleotide shuffling within a sliding window generator
    	d		di-nucleotide shuffling generator
    	w		di-nucleotide shuffling within a sliding window generator
    	g		%GC distribution-based background chooser
    	c		%GC distribution and %GC composition within a sliding window background chooser

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit

To view the help for the individual subcommands, please type:

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
