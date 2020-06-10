============
Introduction
============

The BiasAway software tool is introduced to generate nucleotide
composition-matched DNA sequences. It is available as open source code from
bitbucket.

The tool provides users with six approaches to generate synthetic or genomic
background sequences matching mono- and k-mer composition of user-provided
foreground sequences:

 1) synthetic mononucleotide shuffled sequences
 2) synthetic k-mer shuffled sequences
 3) synthetic mononucleotide shuffled sequences in a sliding window
 4) synthetic k-mer shuffled sequences in a sliding window
 5) genomic mononucleotide distribution matched sequences
 6) genomic mononucleotide distribution within a sliding window matched sequences

The 1st and 2nd approaches shuffle each user-provided sequences independently
by preserving the mononucleotide or k-mer composition, respectively. The 3rd
and 4th approaches apply the same method as for the 1st and 2nd approaches but
within a sliding window along the user-provided sequences. For the 5th and 6th
approaches, the background sequences are selected from a pool of provided
genomic sequences to match the distribution of mononucleotide for each target
sequence. The 6th approach consideres the mean and standard deviation of %GC
computed within the sliding window along the user-provided sequences to match
as closely as possible the distribution for each user-provided sequence.

The approaches based on a sliding window were considered because due to
evolutionary changes such as insertion of repetitive sequences, local
rearrangements, or biochemical missteps, the target sequences may have
sub-regions of distinct nucleotide composition.
