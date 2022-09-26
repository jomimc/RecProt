"""
Tools for studying enzymes that process amino acids.
"""
from matplotlib.lines  import Line2D
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns

protname2letter = {'Alanine':'A', 'Cystine':'C', 'Aspartic acid':'D', 'Glutamic acid':'E', 'Phenylalanine':'F', 'Glycine':'G', 'Histidine':'H', 'Isoleucine':'I', 'Lysine':'K', 'Leucine':'L', 'Methionine':'M', 'Asparagine':'N', 'Proline':'P', 'Glutamine':'Q', 'Arginine':'R', 'Serine':'S', 'Threonine':'T', 'Valine':'V', 'Tryptophan':'W', 'Tyrosine':'Y', 'Aspartate':'D', 'Glutamate':'E', 'Cysteine':'C', 'Glutamyl':'E'}


def load_uniprot_ARS(db='swiss'):
    if db == 'swiss':
        df = pd.read_csv('swissprot_trna_synthetase.tsv', delimiter='\t')
    else:
        df = pd.read_csv('trembl_trna_synthetase.tsv', delimiter='\t')
    df['cognate'] = df['Protein names'].apply(lambda x: protname2letter.get(x.split('-')[0], ''))
    df = df.loc[df.cognate!='']
    df = df.loc[df['Protein families'].notnull()]
    df['fam'] = df['Protein families'].apply(parse_class)
    return df


def load_nonAARS_enzymes():
    df = pd.read_pickle('amino_acid_substrate.pkl')
    df2 = pd.read_excel('tRNA_specificity.xlsx', sheet_name='enzymes')
    ec_key = {ec:c for ec, c in zip(df2['EC number'], df2.Cognate)}
    df['Cognate'] = df.EC.apply(lambda x: ec_key.get(x, ''))
    return df


def load_both_ARS_non_ARS():
    df1 = load_uniprot_ARS()
    df2 = load_nonAARS_enzymes()
    
    L = np.append(df1.Length.values, df2.AA.values)
    T = ['AARS'] * len(df1) + ['non-AARS'] * len(df2)
    C = np.append(df1.cognate.values, df2.Cognate.values)

    df = pd.DataFrame(data={"enzyme_type":T, "Length":L, 'cognate':C})
    df.Length = df.Length.astype(int)
    return df


def parse_class(txt):
    if txt[:3] == 'Cla':
        return txt.split()[0]
    else:
        a = txt.split()[-2]
        if a == '1':
            return "Class-I"
        else:
            return "Class-II"


def print_sig_matched_AA(df):
    AA = np.array(sorted(df.cognate.unique()))
    i0 = df.enzyme_type=='AARS'
    i1 = df.enzyme_type!='AARS'
    for a in AA:
        if a in df.loc[i0, 'cognate'].values and a in df.loc[i1, 'cognate'].values:
            print(a)
            print('AARS:', np.mean(df.loc[(i0)&(df.cognate==a), "Length"]))
            print('non-AARS:', np.mean(df.loc[(i1)&(df.cognate==a), "Length"]))
            print(mannwhitneyu(*[df.loc[(i)&(df.cognate==a), "Length"].values for i in [i0,i1]]))


def plot_length_diff(df):
    fig, ax = plt.subplots()
    sns.boxenplot(x='enzyme_type', y='Length', data=df)
    ax.set_xlabel("Enzyme Type")
    ax.set_ylabel("Sequence Length")
    fig.savefig("si4.pdf", bbox_inches='tight')


def plot_selectivity():
    df = pd.read_excel('tRNA_specificity.xlsx', sheet_name='substrate_selectivity')
    df['log_selectivity'] = np.log10(df.Selectivity)
    sns.catplot(x="is_AARS", y="log_selectivity", data=df)


def get_mean_len(df, AA):
    return np.array([df.loc[(df.cognate==a), 'Length'].mean() for a in AA])


def get_modal_class(classes):
    cuniq = np.unique(classes)
    count = np.array([np.sum(classes==c) for c in cuniq])
    return cuniq[np.argmax(count)]


def fig6(df):
#   df = load_uniprot_ARS(db='swiss')
#   df = df.loc[df.subunit==False]
    AA = np.array(sorted(df.cognate.unique()))
    palette = [sns.color_palette("colorblind")[2], sns.color_palette("colorblind")[3]]
    lbls = ["Class-I", "Class-II"]
    color_key = {l:p for l, p in zip(lbls, palette)}
    fs = 14

    idx = []
    col = []
    for a in AA:
        Class = get_modal_class(df.loc[df.cognate==a, 'fam'])
        lo, hi = df.loc[(df.cognate==a), 'Length'].quantile([0.005, 0.995])
        idx.extend(list(df.loc[(df.cognate==a)&(df.Length>=lo)&(df.Length<=hi)].index))
        col.append(color_key[Class])
    df = df.loc[idx]

    mean_len = get_mean_len(df, AA)
    order = AA[np.argsort(mean_len)[::-1]]

    fig, ax = plt.subplots(figsize=(12,4))
#   sns.boxenplot(x='cognate', y='Length', data=df,
    sns.violinplot(x='cognate', y='Length', data=df, scale='width',
                order=order, showfliers=False, palette=col)
    ax.set_ylim(0, 1300)
    post_cat_edit = ["A", "F", "T", "V", "I", "L", "P"]
    for l in ax.get_xticklabels():
        l.set_fontsize(fs)
        if l.get_text() in post_cat_edit:
            l.set_color(sns.color_palette("dark")[3])
    ax.set_xlabel("Cognate Amino Acid", fontsize=fs)
    ax.set_ylabel("Sequence Length", fontsize=fs)

#   ax.set_yticks()
    for l in ax.get_yticklabels():
        l.set_fontsize(fs)

    handles = [Patch([], [], color=c, ec='k') for c in palette]
    ax.legend(handles, lbls, frameon=False, loc='upper right', fontsize=fs)

    fig.savefig("/home/johnmcbride/LaTEX/RecProtPaper/Figures/fig7.pdf", bbox_inches='tight')


def ligand_pairing():
    df = pd.read_excel('tRNA_specificity.xlsx', sheet_name='substrate_selectivity')
    df = df.loc[df.is_AARS=='T']
    pair_list = [[("Val", "Ile"), ("Ser", "Thr")],
             [("Ala", "Ser"), ("Phe", "Tyr")]]

    size_list = [[[951, 939], [860, 1284]],
                 [[1752, 860], [2244, 848]]]

    xticklabels = [('Smaller by\n-CH3', 'Larger by\n-CH3'),  ('One less\n-OH', 'One more\n-OH')]
    col = [sns.color_palette("colorblind")[2], sns.color_palette("colorblind")[3]]

    fig, ax = plt.subplots(1,2, figsize=(6,2.5))
#   fig, ax = plt.subplots(1,2, figsize=(12,5.5))
    ax1 = [a.twinx() for a in ax]
    fig.subplots_adjust(wspace=0.8)
    for i, pairs in enumerate(pair_list):
        for j, (c, nc) in enumerate(pairs):
            forward = df.loc[(df.Cognate==c)&(df['Non-cognate']==nc), 'Selectivity'].mean()
            backward = df.loc[(df.Cognate==nc)&(df['Non-cognate']==c), 'Selectivity'].mean()
            ax[i].plot([forward, backward], '-o', label=f"Sel, {c}-{nc}", c=col[j])
            ax1[i].plot([0,1], size_list[i][j], ':o', label=f"Size, {c}-{nc}", c=col[j])
        ax[i].set_yscale('log')
        ax[i].set_xlabel("Cognate Ligand")
        ax[i].set_xticks([0,1])
        ax[i].set_xticklabels(xticklabels[i])
        ax[i].set_ylabel("Pairwise Selectivity")
        ax1[i].set_ylabel("Size of ARS complex")

#       ax[i].legend(loc='upper right', bbox_to_anchor=(0.45, 1.25), frameon=False)
#       ax1[i].legend(loc='upper left', bbox_to_anchor=(0.55, 1.25), frameon=False)
        ax[i].set_xlim(-.3, 1.3)

#       ax[i].spines['top'].set_visible(False)
#       ax[i].spines['right'].set_visible(False)
#       ax1[i].spines['top'].set_visible(False)
#       ax1[i].spines['right'].set_visible(False)

    h1 = [Line2D([], [], ls='-', c=c) for c in col]
    h2 = [Line2D([], [], ls=p, c='k') for p in '-:']
    l1 = ["Val-Ile", "Ser-Thr"]
    l2 = ["Ala-Ser", "Phe-Tyr"]
    l3 = ["Sel", "Size"]

    ax[0].legend(h1, l1, bbox_to_anchor=(0.75, 1.25), frameon=False, ncol=1)
    ax1[0].legend(h2, l3, bbox_to_anchor=(1.75, 1.25), frameon=False, ncol=1)
    ax[1].legend(h1, l2, bbox_to_anchor=(0.85, 1.25), frameon=False, ncol=1)

    fig.savefig("fig7.pdf", bbox_inches='tight')
            

def ligand_pairing_2():
    df = pd.read_excel('tRNA_specificity.xlsx', sheet_name='substrate_selectivity')
    df = df.loc[df.is_AARS=='T']
    pair_list = [[("Thr", "Ser"), ("Ile", "Val")],
             [("Phe", "Tyr"), ("Ala", "Ser")]]

    size_list = np.array([1284, 860, 939, 951, 2244, 848, 1752, 860])

    col = [sns.color_palette("colorblind")[2], sns.color_palette("colorblind")[3]]

    fig, ax = plt.subplots(figsize=(6,3.0))
#   fig.subplots_adjust(wspace=0.8)
    X = np.append(np.arange(4), np.arange(5,9))
    X = np.array([0, 1, 2.5, 3.5, 5.5, 6.5, 8, 9])
    width = 0.3
    pair_sel = []
    ax1 = ax.twinx()
    for i, pairs in enumerate(pair_list):
        for j, (c, nc) in enumerate(pairs):
            forward = df.loc[(df.Cognate==c)&(df['Non-cognate']==nc), 'Selectivity'].mean()
            backward = df.loc[(df.Cognate==nc)&(df['Non-cognate']==c), 'Selectivity'].mean()
            pair_sel.extend([forward, backward])
    ax.bar(X-width/2, pair_sel, width, color=col[0], ec='k', lw=1)
    ax1.bar(X+width/2, size_list, width, color=col[1], ec='k', lw=1)
    ax.set_yscale('log')
    ax.set_xticks(X)
    ax.set_xticklabels([x for y in pair_list for z in y for x in z])
    ax.set_ylabel("Pairwise Selectivity")
    ax1.set_ylabel("Size of ARS complex")

    ax.spines['top'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    xpnts = [-.5, 1.75, 4, 0.75]
    xp = np.array([xpnts[0]]*2 + [xpnts[1]]*3 + [xpnts[2]]*2)
    ypnts = [1500, 1600, 1700, 1800]
    yp = np.array([ypnts[0]] + [ypnts[1]]*2 + [ypnts[2]] + [ypnts[1]]*2 + [ypnts[0]])
    ax1.plot(xp, yp, '-k', lw=1.5)
    ax1.text(xpnts[3], ypnts[3], 'Differ by -CH3', fontsize=12)

    dx, dy = 5.5, 700
    ax1.plot(xp + dx, yp + dy, '-k', lw=1.5)
    ax1.text(xpnts[3] + dx, ypnts[3] + dy, 'Differ by -OH', fontsize=12)

    ax.set_xlabel("Cognate Ligand")
    ax.set_ylim(90, 4000000)

    handles = [Patch([], [], color=c, ec='k') for c in col]
    lbls = ["Sel", "Size"]

    ax.legend(handles, lbls, bbox_to_anchor=(0.4, 1.05), frameon=False, ncol=2)

    fig.savefig("/home/johnmcbride/LaTEX/RecProtPaper/Figures/fig7.pdf", bbox_inches='tight')
            

def ligand_pairing_3():
    df = pd.read_excel('tRNA_specificity.xlsx', sheet_name='substrate_selectivity')
    df = df.loc[df.is_AARS=='T']
    pair_list = [[("Thr", "Ser"), ("Ile", "Val")],
             [("Phe", "Tyr"), ("Ala", "Ser")]]

    size_list = -np.array([1284, 860, 939, 951, 2244, 848, 1752, 860])

    col = [sns.color_palette("colorblind")[2], sns.color_palette("colorblind")[3]]

    fig, ax = plt.subplots(2,1,figsize=(5.4,5.4))
    fig.subplots_adjust(hspace=0.45)
    X = np.append(np.arange(4), np.arange(5,9))
    X = np.array([0, 1, 2.5, 3.5, 5.5, 6.5, 8, 9])
    width = 0.3
    pair_sel = []
    for i, pairs in enumerate(pair_list):
        for j, (c, nc) in enumerate(pairs):
            forward = df.loc[(df.Cognate==c)&(df['Non-cognate']==nc), 'Selectivity'].mean()
            backward = df.loc[(df.Cognate==nc)&(df['Non-cognate']==c), 'Selectivity'].mean()
            pair_sel.extend([forward, backward])
    ax[0].bar(X-0*width/2, pair_sel, width, color=col[0], ec='k', lw=1)
    ax[1].bar(X+0*width/2, size_list, width, color=col[1], ec='k', lw=1)
    ax[0].set_yscale('log')
    ax[0].set_xticks(X)
    ax[1].set_xticks(X)
    ax[0].set_xticklabels([x for y in pair_list for z in y for x in z])
    ax[1].set_xticklabels([x for y in pair_list for z in y for x in z])
#   ax[1].set_xticklabels([''] * len(X))
    ax[0].set_ylabel("Pairwise Selectivity")
    ax[1].set_ylabel("Size of ARS complex \n (amino acids)")

    yticks = np.arange(0, -2200, -500)
    ax[1].set_yticks(yticks)
    ax[1].set_yticklabels(-yticks)

    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    xpnts = [-.5, 1.75, 4, 0.75]
    xp = np.array([xpnts[0]]*2 + [xpnts[1]]*3 + [xpnts[2]]*2)
#   ypnts = [1500, 1600, 1700, 1800]
    ypnts = np.array([8000. ** (1.05**(i)) for i in range(4)])
    yp = np.array([ypnts[0]] + [ypnts[1]]*2 + [ypnts[2]] + [ypnts[1]]*2 + [ypnts[0]])
    ax[0].plot(xp, yp, '-k', lw=1.5)
    ax[0].text(xpnts[3], ypnts[3], 'Differ by -CH3', fontsize=12)

#   dx, dy = 5.5, 700
    dx, dy = 5.5, 15.0
    ax[0].plot(xp + dx, yp * dy, '-k', lw=1.5)
    ax[0].text(xpnts[3] + dx, ypnts[3] * dy, 'Differ by -OH', fontsize=12)

    ax[0].set_xlabel("Cognate Ligand", weight='bold')
    ax[0].set_xlim(-1, 9.5)
    ax[1].set_xlim(-1, 9.5)
    ax[0].set_ylim(90, 1000000)
    ax[1].xaxis.set_label_position("top")
    ax[1].xaxis.set_ticks_position("top")

    for i in [0,1]:
        for j in range(0,8,2):
            ax[i].get_xticklabels()[j].set_color("red")

#   handles = [Patch([], [], color=c, ec='k') for c in col]
#   lbls = ["Sel", "Size"]
#   ax[0].legend(handles, lbls, bbox_to_anchor=(0.4, 1.05), frameon=False, ncol=2)

    fig.savefig("/home/johnmcbride/LaTEX/RecProtPaper/Figures/fig7.pdf", bbox_inches='tight')
            

def shear_vs_dist(df):
#   df = pd.read_pickle("/home/johnmcbride/projects/RecProt/ExpData/shear_vs_dist.pkl")
#   fig, ax = plt.subplots(figsize=(12,6))
    fig, ax = plt.subplots(8,1, sharex=True, sharey=True, figsize=(6,12))
    bins = np.arange(-6.05, 0.1, 0.1)
    X = np.arange(-6, 0.1, 0.1)
    col = sns.color_palette('deep')
    for i in range(8):
#       sns.distplot(df.loc[df.dx==(i+1)*5, 'log_shear'], ax=ax[i], bins=)
        hist = np.histogram(df.loc[df.dx==(i+1)*5, 'log_shear'], density=True, bins=bins)[0]
        cum_dist = np.cumsum(hist / hist.sum())
        imax = np.argmin(np.abs(cum_dist-0.5))
        lbl = r"${0:d}<dx\leq{1:d}$ nm".format(i*5, (i+1)*5)
        ax[i].plot(X, hist, label=lbl, c=col[i])
        ax[i].plot([X[imax]]*2, [0, hist[imax]], '-k')
        ax[i].set_xlim(-6, 0)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].legend(loc='upper right', frameon=False)
        if i == 7:
            ax[i].set_xlabel(r"$\log_{10}$ Shear")
#       if i % 2:
        ax[i].set_ylabel("Density")
#   ax.set_xticks(range(8))
#   ax.set_xticklabels([r"${0:d}<dx\leq{1:d}$".format(i*5, (i+1)*5) for i in range(8)])
#   ax.set_xlabel(r"$dx$")
#   ax.set_ylabel(r"$\log_{10}$ Shear")

    fig.savefig("/home/johnmcbride/LaTEX/RecProtPaper/Figures/si3_3.pdf", bbox_inches='tight')


if __name__ == "__main__":


    df = load_both_ARS_non_ARS()
    plot_length_diff(df)
    ligand_pairing_3()


