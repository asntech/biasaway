"""
Modified by Aziz Khan and Anthony Mathelier
"""

from Bio import SeqIO
from Bio.Data import IUPACData
from os.path import splitext
import gzip
import itertools
import seaborn as sns

from sklearn.metrics import mean_absolute_error as mae

IUPAC = list(IUPACData.ambiguous_dna_letters)
IUPAC_DINUC = [''.join(letters) for letters in itertools.product(IUPAC,
                                                                 repeat=2)]


def open_for_parsing(filename):
    """
    Open the file with the given filename for parsing. stdin is used for
    reading from the standard Unix input. Where the file appears to be a
    gzip-compressed file it is opened in the appropriate way to provide access
    to the decompressed text.
    """
    ext = splitext(filename)[-1]
    if ext == ".gz":
        return gzip.open(filename, "rt")  # text flag "t" needed in Python 3
    elif filename == "stdin":
        import sys
        return sys.stdin
    else:
        return open(filename)


def get_seqs(f):
    """
    Retrieve sequences from the input file in fasta format (can be gzipped).
    """
    seqs = []
    fg_gc_list = []
    fg_lengths = []
    dinuc = [0] * len(IUPAC) * len(IUPAC)
    # We cannot use 'with' statement as it would otherwise close stdin
    # automatically
    stream = open_for_parsing(f)
    for record in SeqIO.parse(stream, "fasta"):
        record.seq = record.seq.upper()
        seqs.append(record)
        fg_gc_list.append(GC(record.seq))
        fg_lengths.append(len(record.seq))
        dinuc = [x + y for x, y in zip(dinuc, dinuc_count(record.seq))]
    # Not that we are done, we need to close the opened file but make sure that
    # we do not close stdin
    if f != "stdin":
        stream.close()
    return seqs, fg_gc_list, fg_lengths, dinuc


def GC(seq):
    """
    Calculate G+C content, return the percentage (float between 0 and 100).
    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content. The percentage is calculated against
    the length of the sequence using the IUPAC alphabet, e.g.:
    >>> GC("ACTGN")
    50.0 Note that this will return zero for an empty sequence.

    """
    try:
        gc = sum(map(seq.count, ['G', 'C', 'g', 'c', 'S', 's']))
        length = sum(map(seq.count, ['G', 'C', 'A', 'T', 'S', 'W', 'R', 'Y',
                                     'K', 'M', 'B', 'D', 'H', 'V', 'N', 'g',
                                     'c', 'a', 't', 's', 'w', 'r', 'y', 'k',
                                     'm', 'b', 'd', 'h', 'v', 'n']))
        # fix for py3 issues
        return int(round(gc * 100 / length))
    except ZeroDivisionError:
        return 0


def dinuc_count(seq):
    """
    Calculate dinucleotide composition of a sequence.
    We consider the IUPAC alphabet here.
    """
    return map(seq.count_overlap, IUPAC_DINUC)


def above_threshold(the_list, threshold):
    """
    Returns the percentage of values strictly above the given threshold.
    """
    import numpy as np
    return np.count_nonzero(np.array(the_list) < threshold) / len(the_list)
    # return(iterlen(val for val in the_list if val < threshold) /
    # len(the_list))


def power_div(fg_dist, bg_dist, lambda_="pearson"):
    """
    Compute the power divergence between two distributions.
    Need to test for statistical power first. Given scipy documentation:
        This test is invalid when the observed or expected frequencies in each
        category are too small. A typical rule is that all of the observed and
        expected frequencies should be at least 5.
    """
    import numpy as np
    # Removing 0-values in the expectation
    f_exp = np.asanyarray(fg_dist)
    f_obs = np.asanyarray(bg_dist)[f_exp != 0]
    f_exp = f_exp[f_exp != 0]
    if ((above_threshold(f_exp, 5) > 0.2) or
       (above_threshold(f_obs, 5) > 0.2)):
        return None, -1
    from scipy.stats import power_divergence
    # With the scipy implementation, it is expected that the number of
    # observation and expected counts are the same! Hence, we need to computed
    # the proportion of expected counts in each category to estimate the
    # expected numbers of observations.
    f_exp = sum(f_obs) * (f_exp / sum(f_exp))
    return power_divergence(f_exp=f_exp, f_obs=f_obs, lambda_=lambda_)


def single_value(the_array):
    """
    Assert if the array contains a single value.
    """
    import numpy as np
    return len(np.unique(the_array)) == 1


def QC_info(fg_hist, bg_hist, out_filename):
    """
    Compute QC metrics and print info to out_filename.
    """
    mean_abs_error = mae(fg_hist, bg_hist)
    chi_stat, chi_pval = power_div(fg_hist, bg_hist)
    gof_stat, gof_pval = power_div(fg_hist, bg_hist, "cressie-read")
    gte_stat, gte_pval = power_div(fg_hist, bg_hist, "log-likelihood")
    with open(out_filename, 'w') as stream:
        if sum(fg_hist) < 1000 or sum(bg_hist) < 1000:
            stream.write("QC tests cannot be ")
            stream.write("computed due to a small number of samples ")
            stream.write("(less than 1000).\n")
        elif chi_pval == -1 or gof_pval == -1:
            stream.write("QC tests cannot be ")
            stream.write("computed due to a large number of values ")
            stream.write("with low frequencies (more than ")
            stream.write("20% of values <=5).\n")
        else:
            stream.write("mean_absolute_error\t%f\n" % mean_abs_error)
            stream.write("chi-square(statistic, pvalue)\t(%f, %f)\n" %
                         (chi_stat, chi_pval))
            stream.write("cressie-read goodness_of_fit(statistic, pvalue)")
            stream.write("\t(%f, %f)\n" % (gof_stat, gof_pval))
            stream.write("G-test goodness_of_fit(statistic, pvalue)")
            stream.write("\t(%f, %f)\n" % (gte_stat, gte_pval))


def make_gc_plot(fg_gc, bg_gc, plot_filename):
    """
    Compute the density GC composition plots for background and input.
    """
    from numpy import histogram
    import matplotlib
    matplotlib.use('Agg')
    fg_hist, _ = histogram(fg_gc, bins=101, range=(0.0, 100.))
    bg_hist, _ = histogram(bg_gc, bins=101, range=(0.0, 100.))
    plot_hist = False
    plot_kde = True
    ylab = "density"
    # One cannot compute kde if there is a single value in the array
    if single_value(fg_gc) or single_value(bg_gc):
        plot_hist = True
        plot_kde = False
        ylab = "frequency"
    import matplotlib.pyplot as plt
    plt.figure()
    plot = sns.distplot(fg_gc, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='input')
    plot = sns.distplot(bg_gc, hist=plot_hist, kde=plot_kde,
                        kde_kws={'shade': True, 'linewidth': 3},
                        label='generated')
    plt.legend()
    plot.set(xlabel="%GC", ylabel=ylab)
    plt.savefig("{0}_gc_plot.png".format(plot_filename))
    QC_info(fg_hist, bg_hist, "{0}_gc_plot_stats.txt".format(plot_filename))


def make_len_plot(fg_len, bg_len, plot_filename):
    """
    Compute the density length plot for the input and the matching background
    datasets.
    """
    from numpy import histogram
    from numpy import hstack
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    min_len = min(min(fg_len), min(bg_len))
    max_len = max(max(fg_len), max(bg_len))
    fg_hist, _ = histogram(fg_len, bins=max_len - min_len + 1,
                           range=(min_len, max_len))
    bg_hist, _ = histogram(bg_len, bins=max_len - min_len + 1,
                           range=(min_len, max_len))
    plot_hist = False
    ylab = "density"
    # One cannot compute kde if there is a single value in the array
    if single_value(fg_len) or single_value(bg_len):
        plot_hist = True
        ylab = "frequency"
    plt.figure()
    if plot_hist:
        bins = histogram(hstack((fg_len, bg_len)), bins=100)[1]
        plt.hist(fg_len, bins=bins, align='right', alpha=.5, label='input')
        plt.hist(bg_len, bins=bins, align='right', alpha=.5, label='generated')
    else:
        sns.distplot(fg_len, hist=plot_hist, kde=True, kde_kws={'shade': True,
                                                                'linewidth':
                                                                3},
                     label='input')
        sns.distplot(bg_len, hist=plot_hist, kde=True, kde_kws={'shade': True,
                                                                'linewidth':
                                                                3},
                     label='generated')
    plt.xlabel('length')
    plt.ylabel(ylab)
    plt.legend()
    basename = "{0}_length_plot".format(plot_filename)
    QC_info(fg_hist, bg_hist, "{0}_stats.txt".format(basename))
    plt.savefig("{0}.png".format(basename))


def make_dinuc_dico(dinuc_counts):
    """
    Create a dictionary with all dinucleotide frequencies.
    """
    dico = {}
    for indx, nuc in enumerate(IUPAC):
        min_indx = indx * len(IUPAC)
        max_indx = min_indx + len(IUPAC)
        dico[nuc] = dinuc_counts[min_indx:max_indx]
    return dico


def make_dinuc_plot(fg_dinuc, bg_dinuc, plot_filename):
    """
    Plot the dinucleotide composition of input and background sequences.
    We use the IUPAC alphabet.
    """
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    QC_info(fg_dinuc, bg_dinuc,
            "{0}_dinuc_plot_stats.txt".format(plot_filename))
    fg_total = sum(fg_dinuc)
    fg_dinuc = [val / fg_total for val in fg_dinuc]
    dico = make_dinuc_dico(fg_dinuc)
    fg_df = pd.DataFrame(dico, index=IUPAC)
    bg_total = sum(bg_dinuc)
    bg_dinuc = [val / bg_total for val in bg_dinuc]
    dico = make_dinuc_dico(bg_dinuc)
    bg_df = pd.DataFrame(dico, index=IUPAC)
    mini = min(min(fg_dinuc), min(bg_dinuc))
    maxi = max(max(fg_dinuc), max(bg_dinuc))
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    fig.subplots_adjust(wspace=0.05)
    sns.heatmap(fg_df, cmap="icefire", ax=ax1, cbar=False,
                vmin=mini, vmax=maxi)  # annot=True, fmt='.2f'
    ax1.set(xlabel="first nucleotide", ylabel="second nucleotide",
            title="input")
    fig.colorbar(ax1.collections[0], ax=ax1, location="left",
                 use_gridspec=False, pad=0.2)
    sns.heatmap(bg_df, cmap="icefire", ax=ax2, cbar=False,
                vmin=mini, vmax=maxi)  # annot=True, fmt='.2f"
    ax2.set(xlabel="first nucleotide", title="generated")
    fig.colorbar(ax2.collections[0], ax=ax2, location="right",
                 use_gridspec=False, pad=0.2)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelrotation=0)
    plt.savefig("{0}_dinuc_plot.png".format(plot_filename))


def make_dinuc_acgt_only_dico(dinuc_counts):
    """
    Plot the dinucleotide composition of input and background sequences.
    Only the frequencies for ACGT are considered here.
    """
    acgt = "ACGT"
    dico = {}
    for nuc in acgt:
        dico[nuc] = {}
        for nuc2 in acgt:
            dico[nuc][nuc2] = 0.0
    for indx, dinuc in enumerate(IUPAC_DINUC):
        if(dinuc[0] in "ACGT" and dinuc[1] in "ACGT"):
            dico[dinuc[0]][dinuc[1]] = dinuc_counts[indx]
    final_dico = {}
    for nuc in acgt:
        final_dico[nuc] = []
        for nuc2 in acgt:
            final_dico[nuc].append(dico[nuc][nuc2])
    return final_dico


def make_dinuc_acgt_only_plot(fg_dinuc, bg_dinuc, plot_filename):
    """
    Plot the dinucleotide composition of input and background sequences.
    We only consider A, C, G, T letters.
    """
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    acgt = ['A', 'C', 'G', 'T']
    fg_dico = make_dinuc_acgt_only_dico(fg_dinuc)
    fg_dinuc_count = [item for sublist in fg_dico.values() for item in sublist]
    bg_dico = make_dinuc_acgt_only_dico(bg_dinuc)
    bg_dinuc_count = [item for sublist in bg_dico.values() for item in sublist]
    basename = "{0}_dinuc_acgt_only_plot".format(plot_filename)
    QC_info(fg_dinuc_count, bg_dinuc_count, "{0}_stats.txt".format(basename))
    fg_total = sum(fg_dinuc)
    fg_dinuc = [val / fg_total for val in fg_dinuc]
    fg_dico = make_dinuc_acgt_only_dico(fg_dinuc)
    fg_df = pd.DataFrame(fg_dico, index=acgt)
    bg_total = sum(bg_dinuc)
    bg_dinuc = [val / bg_total for val in bg_dinuc]
    bg_dico = make_dinuc_acgt_only_dico(bg_dinuc)
    bg_df = pd.DataFrame(bg_dico, index=acgt)
    mini = min(min(fg_dinuc), min(bg_dinuc))
    maxi = max(max(fg_dinuc), max(bg_dinuc))
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    fig.subplots_adjust(wspace=0.05)
    sns.heatmap(fg_df, cmap="icefire", ax=ax1, cbar=False,
                vmin=mini, vmax=maxi, annot=True, fmt='.2f')
    ax1.set(xlabel="first nucleotide", ylabel="second nucleotide",
            title="input")
    fig.colorbar(ax1.collections[0], ax=ax1, location="left",
                 use_gridspec=False, pad=0.2)
    sns.heatmap(bg_df, cmap="icefire", ax=ax2, cbar=False, vmin=mini,
                vmax=maxi, annot=True, fmt='.2f')
    ax2.set(xlabel="first nucleotide", title="generated")
    fig.colorbar(ax2.collections[0], ax=ax2, location="right",
                 use_gridspec=False, pad=0.2)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelrotation=0)
    plt.savefig("{0}.png".format(basename))
