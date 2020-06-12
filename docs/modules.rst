=================
BiasAway modules
=================
The BiasAway software tool is introduced to generate nucleotide
composition-matched DNA sequences. It is available as open source code from
bitbucket.

The tool provides users with four approaches to generate synthetic or genomic
background sequences matching mono- and k-mer composition of user-provided
foreground sequences:

.. note:: BiasAway also comes with a Web App available at http://biasaway.uio.no.


K-mer shuffling
==================================

Each user-provided sequence will be shuffled to keep its k-mer composition.
This module can be used for any k, for instance use -k 1 for conserving the
mononucleotide composition of the input sequences.

**Usage:**

.. code-block:: bash

    biasaway k [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway k --help

**Example:**

.. code-block:: bash

    biasaway k -f path/to/FASTA/file/my_fasta_file.fa

It will output the generated sequences on stdout, keeping the dinucleotide
composition of the input sequence by default (k-mer with k=2 is the default).
If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway d -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-k, --kmer","K-mer to be used for shuffling (``default: 2`` for dinucleotide shuffling)"
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"
     "-e, --seed","Seed number to initialize the random number generator for reproducibility (``default: integer from the current time``)"

K-mer shuffling within a sliding window
================================================

For each user-provided sequence, a window will slide along to shuffle the
nucleotides within the window, keeping the local k-mer composition. As such,
the generated sequences will preserve the local k-mer composition of the input
sequences along them.


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

It will output the generated sequences on stdout, keeping the local
dinucleotide composition of the input sequences (k=2 for dinucleotide shuffling
is used as default). If you wish to save the sequences in a specific file, you
can type:

.. code-block:: bash

    biasaway w -f path/to/FASTA/file/my_fasta_file.fa > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-k, --kmer","K-mer to be used for shuffling (``default: 2`` for dinucleotide shuffling)"
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"
     "-w, --winlen", "Window length (``default: 100``)"
     "-s, --step",  "Sliding step (``default: 1``)"
     "-e, --seed","Seed number to initialize the random number generator for reproducibility (``default: integer from the current time``)"

Genomic mononucleotide distribution matched
============================================

Given a set of available background sequences (pre-computed or provided by the
user), each user-provided foreground sequence will be matched to a background
sequence having the same mononucleotide composition.

**Usage:**

.. code-block:: bash

    biasaway g [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway g --help

**Example:**

.. code-block:: bash

    biasaway g -f path/to/FASTA/file/my_fasta_file.fa -b path/to/background.fa -r path/to/bgdirectory

It will output the generated sequences on stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway g -f path/to/FASTA/file/my_fasta_file.fa -b path/to/background.fa -r path/to/bgdirectory > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"
     "-r, --bgdirectory", "Background directory"
     "-b, --background", "Background file in fasta format"
     "-l, --length", "Try to match the length as closely as possible (``not set by default``)"
     "-e, --seed","Seed number to initialize the random number generator for reproducibility (``default: integer from the current time``)"


Genomic mononucleotide distribution within a sliding window matched
===================================================================

Given a set of available background sequences (pre-computed or provided by the
user), each user-provided foreground sequence will be matched to a background
sequence having a close mononucleotide local composition. Specifically,
distribution of %GC composition in a sliding window are computed for foreground
and background sequences; a foreground sequence with a mean m_f and standard
deviation sdev_f of %GC in the sliding window is matched to a background
sequence if its mean %GC m_b is such that:
.. math::
    m_f - N * sdev_f <= m_b <= m_f + N * sdev_f

with *N* equals to 2.6 by default.

**Usage:**

.. code-block:: bash

    biasaway c [options]

.. note:: Please scroll down to see a detailed summary of available **options**.

**Help:**

.. code-block:: bash

    biasaway c --help

**Example:**

.. code-block:: bash

    biasaway c -f path/to/FASTA/file/my_fasta_file.fa -b path/to/background.fa -r path/to/bgdirectory

It will output the generated sequences on stdout. If you wish to save the sequences in a specific file, you can type:

.. code-block:: bash

    biasaway c -f path/to/FASTA/file/my_fasta_file.fa -b path/to/background.fa -r path/to/bgdirectory > path/to/output/FASTA/file/my_fasta_output.fa

**Summary of options**

.. csv-table::
   :header: "Option", "Description"
   :widths: 10, 80

     "-h, --help","To show the help message and exit"
     "-f, --foreground","Foreground file in fasta format."
     "-n, --nfold","How many background sequences per each foreground sequence will be generated (``default: 1``)"
     "-r, --bgdirectory", "Background directory"
     "-b, --background", "Background file in fasta format"
     "-l, --length", "Try to match the length as closely as possible (``not set by default``)"
     "-w, --winlen", "Window length (``default: 100``)"
     "-s, --step", "Sliding step (``default: 1``)"
     "-d, --deviation", "Deviation from the mean (``default: 2.6 for a threshold of mean + 2.6 * stdev``)"
     "-e, --seed","Seed number to initialize the random number generator for reproducibility (``default: integer from the current time``)"
