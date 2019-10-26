=================
BiasAway modules
=================
BiasAway provides three types of plots to visualize intersections of genomic regions and list sets. These are pairwise heatmap of N genomic region sets, classic Venn diagrams of genomic regions and list sets of up to 6-way and UpSet plots.

.. note:: Hello this is a note!


Mono-nucleotide shuffling
=========================

**Usage:**

.. code-block:: bash

    biasaway m [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway m --help

**Example:**

.. code-block:: bash

    biasaway m -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway m -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"

Mono-nucleotide shuffling within a sliding window
==================================================

**Usage:**

.. code-block:: bash

    biasaway f [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway f --help

**Example:**

.. code-block:: bash

    biasaway f -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway f -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"


Di-nucleotide shuffling generator
==================================

**Usage:**

.. code-block:: bash

    biasaway d [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway d --help

**Example:**

.. code-block:: bash

    biasaway d -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway d -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"


Di-nucleotide shuffling within a sliding window
================================================

**Usage:**

.. code-block:: bash

    biasaway w [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway w --help

**Example:**

.. code-block:: bash

    biasaway w -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway w -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"

%GC distribution-based background
==================================

**Usage:**

.. code-block:: bash

    biasaway g [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway g --help

**Example:**

.. code-block:: bash

    biasaway g -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway g -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"

%GC distribution and %GC composition within a sliding window
=============================================================

**Usage:**

.. code-block:: bash

    biasaway c [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway c --help

**Example:**

.. code-block:: bash

    biasaway c -f path/to/FASTA/file/my_fasta_file.fa

This will output the sequence as stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway c -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"



.. note::  The option ``--htype=tribar``.
