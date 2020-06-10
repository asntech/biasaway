# Luis del Peso
# Modified by A. Mathelier
# Vancouver, Jan 2012

# Modified by Aziz Khan and Anthony Mathelier


from __future__ import print_function
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import IUPACData
from biasaway.utils import GC, dinuc_count, IUPAC_DINUC
from ushuffle import shuffle, set_seed


def shuffle_window(ss, km, wl, step):
    bs = ss[:]
    for i in range(0, len(bs)-1, step):
        shuff_seq = shuffle(str.encode(bs[i:(i+wl)]), km).decode()
        bs = bs[0:i] + shuff_seq + bs[i+wl:]
    return(bs)  # returns shuffled sequence


def generate_sequences(seqs, kmer, winlen, step, nfold, random_seed):
    set_seed(random_seed)
    cpt = 1
    bg_gc_list = []
    bg_lengths = []
    dinuc = [0] * len(IUPAC_DINUC)
    for record in seqs:
        seq = record.seq.__str__()
        descr = "Background sequence for {0:s}".format(record.name, cpt)
        for n in range(0, nfold):
            new_sequence = shuffle_window(seq, kmer, winlen, step)
            new_seq = SeqRecord(Seq(new_sequence,
                                    IUPACData.ambiguous_dna_letters),
                                id="background_seq_{0:d}".format(cpt),
                                description=descr)
            print(new_seq.format("fasta"), end='')
            bg_gc_list.append(GC(new_sequence))
            bg_lengths.append(len(new_sequence))
            dinuc = [x + y for x, y in zip(dinuc, dinuc_count(new_sequence))]
            cpt += 1
    return bg_gc_list, bg_lengths, dinuc
