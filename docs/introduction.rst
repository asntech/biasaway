============
Introduction
============

BiasAway provides user with six approaches for generating a background sequence useful to enrichment analyses. These backgrounds derived from mono- and di- nucleotide shuffled sequences, and genomic sequences matched to the GC content of the target data

 1) mononucleotide shuffled target sequence to preserve the mononucleotide composition of the target sequences,
 2) dinucleotide shuffled target sequence to preserve the dinucleotide composition of the target sequences,
 3) genomic sequences matched to the mononucleotide composition of each target sequence to preserve the non-random association of nucleotides,
 4) sliding window of mononucleotide shuffled target sequence,
 5) sliding window of dinucleotide shuffled target sequence,
 6) genomic sequences matched in windows of internal mononucleotide composition for each target sequence.

The latter three backgrounds (BiasAway 4-6) are variants of the former three backgrounds (BiasAway 1-3), in which we utilized a sliding window over the ChIP-Seq sequences to determine a distribution for local regions of
composition. The background sequence set is then generated (mono- or di-nucleotide shuffle) or selected from a pool of genomic sequences (genomic composition match) to match the distribution of window compositions for each target sequence. These latter backgrounds were considered because due to evolutionary changes such as insertion of repetitive sequences, local rearrangements, or biochemical missteps, the target sequences may have sub-regions of distinct nucleotide composition.

