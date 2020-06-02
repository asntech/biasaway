"""
Modified by Aziz Khan on October 29, 2019
Modified by A. Mathelier in March 2020
"""

from Bio import SeqIO
from os.path import splitext
import gzip
import re


def open_for_parsing(filename):
    """
    Open the file with the given filename for parsing. Where the file appears
    to be a gzip-compressed file it is opened in the appropriate way to provide
    access to the decompressed text.
    """
    ext = splitext(filename)[-1]
    if ext == ".gz":
        return gzip.open(filename, "rt")  # text flag "t" needed in Python 3
    else:
        return open(filename)


def GC(seq):
    """
    Calculate G+C content, return the percentage (float between 0 and 100).
    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content. The percentage is calculated against
    the length of the sequence using A,C,G,T,S,W with Ns, e.g.:
    >>> GC("ACTGN")
    50.0 Note that this will return zero for an empty sequence.

    """
    try:
        gc = sum(map(seq.count, ['G', 'C', 'g', 'c', 'S', 's']))
        length = sum(map(seq.count, ['G', 'C', 'A', 'T', 'S', 'W', 'g', 'c',
                                     'a', 't', 's', 'w']))
        # fix for py3 issues
        return int(round(gc * 100 / length))
    except ZeroDivisionError:
        return 0


def dinuc_count(seq):
    """
    Calculate dinucleotide composition of a sequence.
    Only considers A, C, G, T
    """
    dinuc = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG',
             'GC', 'CA', 'CT', 'CG', 'CC']
    return map(seq.count, dinuc)


def get_seqs(f):
    seqs = []
    fg_gc_list = []
    fg_lengths = []
    dinuc = [0] * 16
    with open_for_parsing(f) as stream:
        for record in SeqIO.parse(stream, "fasta"):
            record.seq = record.seq.upper()
            seqs.append(record)
            fg_gc_list.append(GC(record.seq))
            fg_lengths.append(len(record.seq))
            dinuc = [x + y for x, y in zip(dinuc, dinuc_count(record.seq))]
    return seqs, fg_gc_list, fg_lengths, dinuc


def init_compo(length):
    dico = []
    for i in range(1, length):
        dico.append({})
        j = i - 1
        dico[j]["AA"] = 0.0
        dico[j]["AC"] = 0.0
        dico[j]["AT"] = 0.0
        dico[j]["AG"] = 0.0
        dico[j]["CA"] = 0.0
        dico[j]["CC"] = 0.0
        dico[j]["CT"] = 0.0
        dico[j]["CG"] = 0.0
        dico[j]["GA"] = 0.0
        dico[j]["GC"] = 0.0
        dico[j]["GG"] = 0.0
        dico[j]["GT"] = 0.0
        dico[j]["TA"] = 0.0
        dico[j]["TC"] = 0.0
        dico[j]["TG"] = 0.0
        dico[j]["TT"] = 0.0
    return dico


def length(record):
    return len(record.seq)


def all_pos_dinuc(compo, max_length):
    distrib = {}
    composition = {}
    for key in compo[0].keys():
        cpt = 0.0
        for i in range(1, max_length):
            cpt += compo[i - 1][key]
        composition[key] = cpt
    for k in composition.keys():
        cpt = 4.0  # WARNING: WE DO NOT TAKE ANY 'N' INTO ACCOUNT
        first = k[0]
        for key in composition.keys():
            if key[0] == first:
                cpt += composition[key]
        distrib[k] = (composition[k] + 1.0) / cpt
    return distrib


def pos_by_pos_dinuc(compo, seq_length):
    distrib = init_compo(seq_length)
    for j in range(1, seq_length):
        for k in compo[j - 1].keys():
            if re.search("N", k):
                distrib[j - 1][k] = 0.0
            else:
                cpt = 4.0
                first = k[0]
                for key in compo[j - 1].keys():
                    if key[0] == first:
                        cpt += compo[j - 1][key]
                distrib[j - 1][k] = (compo[j - 1][k] + 1.0) / cpt
    return distrib


def compute_dinuc_distrib(seqs, b=False):
    max_length = max(map(length, seqs))
    compo = init_compo(max_length)
    for i in range(0, len(seqs)):
        seq_length = len(seqs[i].seq)
        for j in range(1, seq_length):
            if seqs[i].seq[j - 1] != 'N' and seqs[i].seq[j] != 'N':
                compo[j - 1]["%s" % (seqs[i].seq[(j - 1):(j + 1)])] += 1.0

    if b:  # Dinucleotide distrib over all positions
        return all_pos_dinuc(compo, max_length)
    else:  # Dinucleotide distrib position by position
        return pos_by_pos_dinuc(compo, seq_length)


def print_dinuc_distrib(dinuc, output):
    import sys
    stream = sys.stdout
    if output:
        stream = open(output, "w")
    for j in range(0, len(dinuc)):
        stream.write("%f, %f, %f, %f, %f, %f, %f, %f, %f," % (dinuc[j]["AA"],
                                                              dinuc[j]["AC"],
                                                              dinuc[j]["AG"],
                                                              dinuc[j]["AT"],
                                                              dinuc[j]["AN"],
                                                              dinuc[j]["CA"],
                                                              dinuc[j]["CC"],
                                                              dinuc[j]["CG"],
                                                              dinuc[j]["CT"]))
        stream.write(" %f, %f, %f, %f, %f, %f, %f, %f, %f," % (dinuc[j]["CN"],
                                                               dinuc[j]["GA"],
                                                               dinuc[j]["GC"],
                                                               dinuc[j]["GG"],
                                                               dinuc[j]["GT"],
                                                               dinuc[j]["GN"],
                                                               dinuc[j]["TA"],
                                                               dinuc[j]["TC"],
                                                               dinuc[j]["TG"]))
        stream.write(" %f, %f, %f, %f, %f, %f, %f\n" % (dinuc[j]["TT"],
                                                        dinuc[j]["TN"],
                                                        dinuc[j]["NA"],
                                                        dinuc[j]["NC"],
                                                        dinuc[j]["NG"],
                                                        dinuc[j]["NT"],
                                                        dinuc[j]["NN"]))
    stream.close()


def compute_nt_distrib(seqs):
    cpt = 4.0
    distrib = {}
    for l in "ACGT":
        distrib[l] = 1.0
    for seq in seqs:
        for letter in seq:
            if letter != 'N':
                distrib[letter] += 1.0
                cpt += 1.0
    for letter in "ACGT":
        distrib[letter] /= cpt
    return distrib


def split_seq(seq):
    return re.split('([!ACGT]+)', seq)


def single_value(the_array):
    import numpy as np
    return len(np.unique(the_array)) == 1


def test_dir(plot_filename):
    """ Test if the directory exists or can be created. """
    import os
    import errno
    try:
        os.makedirs(plot_filename)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(plot_filename):
            pass
        else:
            raise


def make_gc_plot(fg_gc, bg_gc, plot_filename):
    """
    Compute the density GC composition plots for background and input.
    """
    test_dir(plot_filename)
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt
    plot_hist = False
    plot_kde = True
    ylab = "density"
    # One cannot compute kde if there is a single value in the array
    if single_value(fg_gc) or single_value(bg_gc):
        plot_hist = True
        plot_kde = False
        ylab = "frequency"
    plt.figure()
    plot = sns.distplot(fg_gc, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='input')
    plot = sns.distplot(bg_gc, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='generated')
    plt.legend()
    plot.set(xlabel="%GC", ylabel=ylab)
    plt.savefig("{0}/{0}_gc_plot.png".format(plot_filename))


def make_len_plot(fg_len, bg_len, plot_filename):
    """
    Compute the density length plot for the background, the input and the
    matching background datasets.
    """
    test_dir(plot_filename)
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt
    plot_hist = False
    plot_kde = True
    ylab = "density"
    # One cannot compute kde if there is a single value in the array
    if single_value(fg_len) or single_value(bg_len):
        plot_hist = True
        plot_kde = False
        ylab = "frequency"
    plt.figure()
    plot = sns.distplot(fg_len, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='input')
    plot = sns.distplot(bg_len, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='generated')
    plt.legend()
    plot.set(xlabel="length", ylabel=ylab)
    plot.get_figure().savefig("{0}/{0}_length_plot.png".format(plot_filename))


def make_dinuc_plot(fg_dinuc, bg_dinuc, plot_filename):
    """
    Plot the dinucleotide composition of input and background sequences.
    """
    test_dir(plot_filename)
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt
    nuc = ['A', 'C', 'G', 'T']
    fg_total = sum(fg_dinuc)
    fg_dinuc = [val / fg_total for val in fg_dinuc]
    fg_df = pd.DataFrame({'A': fg_dinuc[0:4], 'C': fg_dinuc[4:8],
                          'G': fg_dinuc[8:12], 'T': fg_dinuc[12:16]},
                         index=nuc)
    bg_total = sum(bg_dinuc)
    bg_dinuc = [val / bg_total for val in bg_dinuc]
    bg_df = pd.DataFrame({'A': bg_dinuc[0:4], 'C': bg_dinuc[4:8],
                          'G': bg_dinuc[8:12], 'T': bg_dinuc[12:16]},
                         index=nuc)
    mini = min(min(fg_dinuc), min(bg_dinuc))
    maxi = max(max(fg_dinuc), max(bg_dinuc))
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    fig.subplots_adjust(wspace=0.05)
    sns.heatmap(fg_df, cmap="icefire", ax=ax1, cbar=False, annot=True,
                fmt='.2f', vmin=mini, vmax=maxi)
    ax1.set(xlabel="first nucleotide", ylabel="second nucleotide",
            title="input")
    fig.colorbar(ax1.collections[0], ax=ax1, location="left",
                 use_gridspec=False, pad=0.2)
    sns.heatmap(bg_df, cmap="icefire", ax=ax2, cbar=False, annot=True,
                fmt='.2f', vmin=mini, vmax=maxi)
    ax2.set(xlabel="first nucleotide", title="generated")
    fig.colorbar(ax2.collections[0], ax=ax2, location="right",
                 use_gridspec=False, pad=0.2)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelrotation=0)
    plt.savefig("{0}/{0}_dinuc_plot.png".format(plot_filename))
