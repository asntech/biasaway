"""
Generation of sequences using mono-nucleotide of input sequences.

Module allowing the generation of sequences by using a mono-nucleotide.
shuffling of the given sequences.

Modified by Aziz Khan on October 29, 2019
Modified by A. Mathelier in March 2020

"""

from __future__ import print_function
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from biasaway.utils import GC, dinuc_count
import random


def generate_sequences(seqs, nfold):
    """
    Generate sequences by shuffling input (mononucleotide).

    return tuple containing %GC compo and length distrib of output.

    """

    cpt = 1
    bg_gc_list = []
    bg_lengths = []
    dinuc = [0] * 16
    for record in seqs:
        seq = record.seq.__str__()
        descr = "Background sequence for {0:s}".format(record.name)
        for _ in range(0, nfold):
            new_sequence = ''.join(random.sample(seq, len(seq)))
            new_seq = SeqRecord(Seq(new_sequence, generic_dna),
                                id="background_seq_{0:d}".format(cpt),
                                description=descr)
            print(new_seq.format("fasta"), end='')
            bg_gc_list.append(GC(new_sequence))
            bg_lengths.append(len(new_sequence))
            dinuc = [x + y for x, y in zip(dinuc, dinuc_count(new_sequence))]
            cpt += 1
    return bg_gc_list, bg_lengths, dinuc
