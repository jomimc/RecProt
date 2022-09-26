import os
import pickle

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
from palettable.scientific.diverging import Berlin_20, Roma_20, Berlin_10, Vik_20
import pandas as pd
import seaborn as sns



def gap_vs_kappa_vs_E(df, xlbl='logK', ylbl='dE', zlbl='gap', nbin=10):
    X = df[xlbl].values
    Y = df[ylbl].values
    Z = df[zlbl].values

#   xbins = np.linspace(X.min(), X.max(), nbin + 1)[1:]
#   ybins = np.linspace(Y.min(), Y.max(), nbin + 1)[1:]

    xbins = np.quantile(X, np.linspace(0, 1, nbin+1))[1:]
    ybins = np.quantile(Y, np.linspace(0, 1, nbin+1))[1:]

    print(xbins[0], xbins[-1])
    print(ybins[0], ybins[-1])

    xi = np.digitize(X, xbins)
    yi = np.digitize(Y, ybins)

    Zi = np.zeros((xbins.size, ybins.size), float)
    for i in range(xbins.size):
        for j in range(ybins.size):
            Zi[i,j] = Z[(xi==i)&(yi==j)].mean()
#           idx = np.where((xi==i)&(yi==j))[0]
#           Zi[i,j] = Z[idx].mean() if len(idx) > 10**5 else np.nan

    plot_heatmap(Zi, xbins, ybins, xlbl, ylbl)


def hist2d(df, xlbl='logK', ylbl='dE', nbin=10):
    X = df[xlbl].values
    Y = df[ylbl].values

#   xbins = np.linspace(X.min(), X.max(), nbin + 1)[1:]
#   ybins = np.linspace(Y.min(), Y.max(), nbin + 1)[1:]

    xbins = np.quantile(X, np.linspace(0, 1, nbin+1))[1:]
    ybins = np.quantile(Y, np.linspace(0, 1, nbin+1))[1:]

    xi = np.digitize(X, xbins)
    yi = np.digitize(Y, ybins)

    Zi = np.zeros((xbins.size, ybins.size), float)
    for i in range(xbins.size):
        for j in range(ybins.size):
            Zi[i,j] = np.sum((xi==i)&(yi==j))

    plot_heatmap(Zi, xbins, ybins, xlbl, ylbl)


def plot_heatmap(Zi, xbins, ybins, xlbl, ylbl, fig='', ax=''):
    if isinstance(ax, str):
        fig, ax = plt.subplots()
    im = ax.imshow(Zi.T)
    ax.invert_yaxis()
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_xticks(range(xbins.size)[::4])
    ax.set_yticks(range(ybins.size)[::4])
    ax.set_xticklabels(np.round(xbins, 0).astype(int)[::4], rotation=90)
    ax.set_yticklabels(np.round(ybins, 0).astype(int)[::4])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_combo(bins, combo, totals, recog, typ='dens', istart=0, istep=1, save=False, itst=0):
    fig, ax = plt.subplots()
    xlbl, ylbl = combo
    if typ in ['prob', 'tot']:
        cs = f"{xlbl}_{ylbl}_dens"
    else:
        cs = f"{xlbl}_{ylbl}_{typ}"
    cs0 = f"{xlbl}_{ylbl}_{typ}"
    if typ == "dens":
        R = np.nansum(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
        plot_heatmap(R, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax)
    elif typ == "tot":
        T = np.nansum(np.array([x for x in totals[cs][istart::istep] if len(x)]), axis=0)
        plot_heatmap(T, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax)
    elif typ == "prob":
        R = np.nansum(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
        T = np.nansum(np.array([x for x in totals[cs][istart::istep] if len(x)]), axis=0)
        plot_heatmap(R / T, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax)
    elif typ == 'Z':
        R = np.nanmean(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
        plot_heatmap(R, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax)

    if save:
        for d in ['All', f'Test{itst:02d}', f"{cs}"]:
            fig.savefig(f"../N13/Figures/{d}/e{istart}_{itst:02d}_{cs0}.pdf", bbox_inches='tight')
            fig.savefig(f"../N13/Figures/{d}/e{istart}_{itst:02d}_{cs0}.png", bbox_inches='tight')



def plot_best_df(df, idx, cmap=Vik_20.mpl_colormap, vmin=1, vmax=2, fig='', ax='', bidx=[1,6,2], ec='', cbar=True, side=[0,1], lig=True):
    if isinstance(fig, str):
        fig, ax = plt.subplots(1,2)
    if cbar:
        for a in ax: 
            sc = a.imshow(np.array([[0.1,0.2],[0.1, 0.2]]), cmap=cmap, vmin=vmin, vmax=vmax)
            a.invert_yaxis()
            sc.set_visible(False)

    col1 = df.loc[idx, "pseq"].mean(axis=0)[:-3]
    col2 = df.loc[idx, "pseq"].mean(axis=0)[-3:]
    for i in side:

        xy = df.loc[idx, f"cprot{i+1}"].mean(axis=0)
        if not isinstance(ec, str):
            ec0 = ec.copy()
            ec0[:,1] = ec0[:,1] - ec0[bidx[1],1] + xy[bidx[1],1]

        for j, (x, c) in enumerate(zip(xy, col1)):
            if not isinstance(ec, str):
#               ax[i].plot([ec0[j,0]], [ec0[j,1]], 'ok', fillstyle='none', ms=6)
#               ax[i].plot([ec0[j,0], x[0]], [ec0[j,1], x[1]], '-k', alpha=0.5)
                ax[i].quiver([ec0[j,0]], [ec0[j,1]], [x[0] - ec0[j,0]], [x[1] - ec0[j,1]], angles='xy', scale_units='xy', scale=1)

            c = (c - vmin) / (vmax - vmin)
            ax[i].plot([x[0]], [x[1]], 'o', c=cmap(c), fillstyle='left', ms=10)

        for x, c in zip(xy[bidx], col2):
            c = (c - vmin) / (vmax - vmin)
            ax[i].plot([x[0]], [x[1]], 'o', c=cmap(c), fillstyle='right', ms=10)

        if lig:
            clig = df.loc[idx, f"clig{i+1}"].iloc[0]
            ax[i].plot(*clig.T, 'sk', ms=12, fillstyle='none', mew=2)

    if cbar:
        fig.colorbar(sc, ax=ax[1], fraction=0.046, pad=0.04)


def plot_rec_space(hist, i0=1, cat="meanK_meanE"):
    fig = plt.figure(figsize=(14,8))
    gs = GridSpec(7,5, width_ratios=[1,0.3,1,0.3,1])
    ax = [fig.add_subplot(gs[1:5,0])] + \
         [fig.add_subplot(gs[i:i+3, j+2]) for i in [0,4] for j in [0,2]]
#   fig, ax = plt.subplots(1,3)
#   sets = ['all', f'bind{i0}', f'rec{i0}']
    sets = sorted(list(hist.keys()))
    xlbl, ylbl = cat.split('_')
    xbins, ybins = hist[sets[0]][cat]['bins']
    Z1, Z2, Z3, Z4, Z5 = [hist[s][cat]['hist'] for s in sets]
    Z1[Z1==0] = np.nan

    Z2 = Z2 / Z1
    Z2[Z2==0] = np.nan
    Z2[np.isfinite(Z2)==False] = np.nan
    
    Z3 = Z3 / Z1
    Z3[Z3==0] = np.nan
    Z3[np.isfinite(Z3)==False] = np.nan

    plot_heatmap(np.log10(Z1), xbins, ybins, xlbl, ylbl, fig, ax[0])
    plot_heatmap(Z2, xbins, ybins, xlbl, ylbl, fig, ax[1])
    plot_heatmap(Z3, xbins, ybins, xlbl, ylbl, fig, ax[2])


    Z4 = Z4 / Z1
    Z4[Z4==0] = np.nan
    Z4[np.isfinite(Z4)==False] = np.nan
    
    Z5 = Z5 / Z1
    Z5[Z5==0] = np.nan
    Z5[np.isfinite(Z5)==False] = np.nan

    plot_heatmap(Z4, xbins, ybins, xlbl, ylbl, fig, ax[3])
    plot_heatmap(Z5, xbins, ybins, xlbl, ylbl, fig, ax[4])


def plot_covariance(df, idx, cmap=Vik_20.mpl_colormap, fig='', ax=''):
    if isinstance(fig, str):
        fig, ax = plt.subplots()

    seq = np.array(list(df.loc[idx, 'pseq']))
    cov = np.cov(seq.T)
    np.fill_diagonal(cov, 0)
    vmax = np.max(np.abs(cov))
    vmin = -vmax

    im = ax.imshow(cov, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.invert_yaxis()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    

def plot_recog_example(df, idx, ec='', bidx=[1,6,2]):
    print(f"Total seq: {np.sum(idx)}")
    fig, ax = plt.subplots(1,3)
    plot_best_df(df, idx, fig=fig, ax=ax[:2], ec=ec, bidx=bidx)
    plot_covariance(df, idx, fig=fig, ax=ax[2])





