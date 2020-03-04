"""
Shuffle input sequences within a sliding window, keeping mononuc compo.

Written by Luis del Peso
Modified by A. Mathelier

Modified by Aziz Khan on October 29, 2019
Modified by A. Mathelier in March 2020

"""

from __future__ import print_function
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from biasaway.utils import GC, dinuc_count
import random


def shuffle_window(ss, wl, step):
    """
    Shuffle sequence in the given window.

    Return the shuffled sequence.

    """

    bs = ss[:]
    for i in range(0, len(bs) - 1, step):
        substring = bs[i:(i + wl)]
        bs = ''.join([bs[0:i],
                      ''.join(random.sample(substring, len(substring))),
                      bs[i + wl:]])
    return bs


def generate_sequences(seqs, winlen, step, nfold):
    """
    Shuffle sequences within a sliding window, keeping mononuc compo.

    Return %GC and length distribution of output sequences.

    """

    cpt = 1
    bg_gc_list = []
    bg_lengths = []
    dinuc = [0] * 16
    for record in seqs:
        sequence = record.seq.__str__()
        descr = "Background sequence for {0:s}".format(record.name, cpt)
        for _ in range(0, nfold):
            new_sequence = shuffle_window(sequence, winlen, step)
            new_seq = SeqRecord(Seq(new_sequence, generic_dna),
                                id="background_seq_{0:d}".format(cpt),
                                description=descr)
            print(new_seq.format("fasta"), end='')
            bg_gc_list.append(GC(new_sequence))
            bg_lengths.append(len(new_sequence))
            dinuc = [x + y for x, y in zip(dinuc, dinuc_count(new_sequence))]
            cpt += 1
    return bg_gc_list, bg_lengths, dinuc
