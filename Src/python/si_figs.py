"""
Code for the figures in the supplement.
Will not run without the correct data, or without correct paths to data.
"""
import os
from pathlib import Path
import pickle
import sys

import matplotlib as mpl
from matplotlib import rc
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from palettable.colorbrewer.qualitative import Paired_12
from palettable.colorbrewer.sequential import Purples_9, YlGn_9, YlOrBr_9, YlGn_6
from palettable.scientific.sequential import Acton_20, Oslo_20, Nuuk_20, LaJolla_20, LaPaz_20, Acton_6
from palettable.scientific.diverging import Berlin_20, Roma_20, Berlin_11, Vik_20, Vik_6
import pandas as pd
from PIL import Image
from scipy.spatial.distance import cdist
from scipy.stats import linregress, pearsonr
import seaborn as sns

import main_figs
import local_utils as utils

import rec_utils
import get_entropy as GE


#ATH_FIG = Path("/home/jmcbride/RecProt/Figures")
PATH_FIG = Path("/data/Jmcbride/RecProt/Figures/Paper")
PATH_DATA = Path("/data/Jmcbride/RecProt/Figures/Paper/Data")
PATH_BASE = Path("/molahome/jmcbride/RecProt")



########################################################
### SI FIG 1
###

### Needs the results of the G-MeCh model (df2)
def fig1(df2, ec, ID='w2h1', cmap=Vik_6, ft=10):
    fig = plt.figure(figsize=(12,9))
    gs = GridSpec(3,11, height_ratios=[1, 0.01, 0.5])#, 0.01, 0.5])
    fig.subplots_adjust(wspace=1.2, hspace=0)
    ax = [fig.add_subplot(gs[2,i*3+i:(i+1)*3+i]) for i in [0,1,2]] + \
         [fig.add_subplot(gs[0,i*3+i:(i+1)*3+i]) for i in [0,1,2]]#+ \
#        [fig.add_subplot(gs[4,i*3+i:(i+1)*3+i]) for i in [0,1,2]]

    df = df2.loc[(df2.ID==ID)&(df2.ang2==0.785398)]

    idx_list = [df.index, (df.gap>3.60), (df.gap>3.800)]
    for i, idx in enumerate(idx_list):
        main_figs.plot_best_df(df, idx, fig=fig, ax=[0,ax[i]], cbar=False, side=[1], lig=False, cmap=cmap.mpl_colormap)
        main_figs.plot_covariance(df, idx, fig=fig, ax=ax[3+i], ft=ft)

    for a in ax[:3]:# + ax[6:]:
        for direction in ['left', 'right', 'top', 'bottom']:
            a.spines[direction].set_visible(False)
        a.set_xticks([])
        a.set_yticks([])
        a.set_xlim(-2.3, 2.3)
#       a.set_ylim(-0.9, 1.5)
        a.set_ylim(-2.3, 2.3)

    xt = np.array([0, 4, 9, 13, 16]) - 0.5
    xl = ["Row 1", "Row 2", "Row 3", "Bind"]

    for a in ax[3:6]:
        a.set_xticks(xt)
        a.set_yticks(xt)
        lx = a.set_xticklabels(xl + [""], rotation=0, horizontalalignment='left', fontsize=9)
        ly = a.set_yticklabels([""] + xl, rotation=0, verticalalignment='bottom')

    norm = mpl.colors.Normalize(vmin=1, vmax=2)
    sc = plt.cm.ScalarMappable(cmap=cmap.mpl_colormap, norm=norm)
    sc.set_array([])
    cbar = fig.colorbar(sc, ax=ax[2], fraction=0.046, pad=0.04)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel("Mean sequence", rotation=270, fontsize=ft)

    fs = 16
    ax[0].text( -0.10, 1.03, 'B', transform=ax[0].transAxes, fontsize=fs)
    ax[3].text( -0.30, 1.03, 'A', transform=ax[3].transAxes, fontsize=fs)


    fig.savefig(PATH_FIG.joinpath(f"si1.pdf"), bbox_inches='tight')


########################################################
### SI FIG 1b 
###


### Examples of deformed proteins
### Needs the results of the G-MeCh model (df2)
### "ec" is the equilibrium configuration
def fig1b(df2, ec, ID='w2h1', cmap=Vik_6, ft=10):
    fig, ax = plt.subplots(4,1,figsize=(3,10))

    df = df2.loc[(df2.ID==ID)&(df2.ang2==0.785398)]
    idx_conf = [5058530, 5059625, 5061726, 5069504]
    for i, j in enumerate(idx_conf):
        print(df.loc[j, 'gap'])
        plot_displacement(ec, *df.loc[j, ['cprot2', 'pseq']], ax[i])

    fs = 16
    for a, lbl in zip(ax, 'ABCD'):
        for direction in ['left', 'right', 'top', 'bottom']:
            a.spines[direction].set_visible(False)
        a.set_xticks([])
        a.set_yticks([])
        a.set_xlim(-2.3, 2.3)
        a.set_ylim(-2.3, 2.3)
        a.text( -0.10, 1.03, lbl, transform=a.transAxes, fontsize=fs)
        
    fig.savefig(PATH_FIG.joinpath(f"si1b.pdf"), bbox_inches='tight')


def transform_vector(x0, y0, dx, dy, scale1=2.0, scale2=0.2):
    theta1 = np.arctan2(dx, dy)
    if theta1 < 0:
        theta2 = theta1 - np.pi / 2
    else:
        theta2 = theta1 + np.pi / 2

    x1 = x0 + np.sin(theta2) * scale2 - dx * scale1 / 2
    y1 = y0 + np.cos(theta2) * scale2 - dy * scale1 / 2

    dx = dx * scale1
    dy = dy * scale1

    return x1, y1, dx, dy


def plot_displacement(ec0, xy, seq, ax='', ms=6):
    if isinstance(ax, str):
        fig, ax = plt.subplots()
#   ax.quiver([ec0[j,0]], [ec0[j,1]], [x[0] - ec0[j,0]], [x[1] - ec0[j,1]], angles='xy', scale_units='xy', scale=1)
    col = [Vik_20.mpl_colors[[3 if s == 1 else 18][0]] for s in seq]
    for j, (x, c) in enumerate(zip(xy, col)):
        ax.plot([ec0[j,0]], [ec0[j,1]], 'o', c=c, mec='grey', mew=1, ms=ms)
        ax.plot([x[0]], [x[1]], 'ok', fillstyle='none', ms=ms)
#       ax.quiver([ec0[j,0]], [ec0[j,1]], [x[0] - ec0[j,0]], [x[1] - ec0[j,1]], angles='xy', scale_units='xy', scale=1)
        x1, y1, dx, dy = transform_vector(ec0[j,0], ec0[j,1], x[0] - ec0[j,0], x[1] - ec0[j,1])
        ax.quiver(x1, y1, dx, dy, angles='xy', scale_units='xy', scale=1)




########################################################
### FIG 5
###

### Needs the results of evolve_function.py
def fig2(out):
    fig, ax = plt.subplots(2,3, figsize=(6,6))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    id_list = ['W1H1', 'W2H1', 'W2H2-', 'W2H2', 'W2H3-', 'W3H1', 'W2H3']
    col = np.array(sns.color_palette())[[0,1,3,2,4]]
    size = np.array([7, 13, 16, 18, 18, 19, 22])
    angles = 2 * (90 - np.arange(30, 95, 5))
    lbls1 = ["Easy Task", "Medium Task", "Hard Task"]
#   for j, (frac, trav, neigh) in enumerate(out):
    df = pd.DataFrame(columns=['phi', 'ID', 'theta', 'frac', 'trav', 'neigh', 'nc'])
    for i in range(3):
        Y1max = out[i][-1][1][0][1:-1]
        Y2max = out[i][-1][2][0][1:-1]
        for j, (frac, trav, neigh, nc) in enumerate(out[i]):
            for k in range(1, len(angles)-1):
                df.loc[len(df)] = [(i+1)*10, id_list[j], angles[k], frac[0][k], trav[0][k], neigh[0][k], nc[0][k]]
#   return df
#           ax[0,i].plot(angles[1:-1], trav[0][1:-1], '-o', label=id_list[j], c=col[j])
#           ax[1,i].plot(angles[1:-1], neigh[0][1:-1], '-o', label=id_list[j], c=col[j])
            if j <=3:
                ax[0,i].plot(angles[1:-1], trav[0][1:-1] / Y1max, '-o', label=id_list[j], c=col[j])
                ax[1,i].plot(angles[1:-1], neigh[0][1:-1] / Y2max, '-o', label=id_list[j], c=col[j])
#               print(i, j, [np.nanmean(x) for x in [trav[0][1:-1] / Y1max, neigh[0][1:-1] / Y2max]])
        ax[0,i].set_title(lbls1[i])
        ax[1,i].set_title(lbls1[i])


    for i in [0,1]:
        for j in [0,1,2]:
            ax[i,j].set_xlabel(r"$\theta_C$")
            ax[i,j].set_xticks(np.arange(20, 110, 20))
#           ax[i,j].set_ylim(-0.03, 0.95)
#           ax[i,j].set_ylim(-0.03, 0.95)
    ax[0,0].set_ylabel("Relative Evolvability")
    ax[1,0].set_ylabel("Relative Robustness")


    for a in ax.ravel():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_ylim(0, 1.1)
    
    ax[0,0].legend(loc='upper right', bbox_to_anchor=(4.0, 1.35), frameon=False, ncol=5)

    fs = 16
    for i, b in enumerate('AB'):
        ax[i,0].text( -0.45, 1.05, b, transform=ax[i,0].transAxes, fontsize=fs)

    fig.savefig(PATH_FIG.joinpath(f"si2.pdf"), bbox_inches='tight')



#ef fig3():
#   df = pd.read_pickle("/mnt/bigram2/jmcbride/ProteinDB/PDB/SimSeq/shear_vs_dist.pkl")
#   fig, ax = plt.subplots(figsize=(12,6))
#   sns.boxplot(x='dx', y='log_shear', data=df.loc[df.dx<=40], showfliers=False)
#   ax.set_xticks(range(8))
#   ax.set_xticklabels([r"${0:d}<dx\leq{1:d}$".format(i*5, (i+1)*5) for i in range(8)])
#   ax.set_xlabel(r"$dx$")
#   ax.set_ylabel(r"$\log_{10}$ Shear")

#   fig.savefig(PATH_FIG.joinpath(f"si3.pdf"), bbox_inches='tight')


### Needs the results of evolve_function.py
def fig5(out):
    fig, ax = plt.subplots(2,3, figsize=(6,6))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    id_list = ['W1H1', 'W2H1', 'W2H2-', 'W2H2', 'W2H3-', 'W3H1', 'W2H3']
    col = np.array(sns.color_palette())[[0,1,3,2,4]]
    size = np.array([7, 13, 16, 18, 18, 19, 22])
    angles = 2 * (90 - np.arange(30, 95, 5))
    lbls1 = ["Easy Task", "Medium Task", "Hard Task"]
#   for j, (frac, trav, neigh) in enumerate(out):
    for i, i0 in enumerate([0, 1, 2]):
        Y1max = out[i0][-1][1][0][1:-1]
        Y2max = out[i0][-1][2][0][1:-1]
        for j, (frac, trav, neigh, nc) in enumerate(out[i0]):
            ax[0,i].plot(angles[1:-1], trav[0][1:-1], '-o', label=id_list[j], c=col[j])
            ax[1,i].plot(angles[1:-1], neigh[0][1:-1], '-o', label=id_list[j], c=col[j])
#           if j <=3:
#               ax[0,i].plot(angles[1:-1], trav[0][1:-1] / Y1max, '-o', label=id_list[j], c=col[j])
#               ax[1,i].plot(angles[1:-1], neigh[0][1:-1] / Y2max, '-o', label=id_list[j], c=col[j])
#               print(i, j, [np.nanmean(x) for x in [trav[0][1:-1] / Y1max, neigh[0][1:-1] / Y2max]])
        ax[0,i].set_title(lbls1[i])
        ax[1,i].set_title(lbls1[i])


    for i in [0,1]:
        for j in [0,1,2]:
            ax[i,j].set_xlabel(r"$\theta_C$")
            ax[i,j].set_xticks(np.arange(20, 110, 20))
            ax[0,j].set_ylim(-0.03, 0.95)
            ax[1,j].set_ylim(-0.03, 0.75)
    ax[0,0].set_ylabel("Evolvability")
    ax[1,0].set_ylabel("Robustness")


    for a in ax.ravel():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    
    ax[0,0].legend(loc='upper right', bbox_to_anchor=(4.0, 1.35), frameon=False, ncol=5)

    fs = 16
    for i, b in enumerate('AB'):
        ax[i,0].text( -0.45, 1.05, b, transform=ax[i,0].transAxes, fontsize=fs)

    fig.savefig(PATH_FIG.joinpath(f"si5.pdf"), bbox_inches='tight')


def bootstrap_cov(ec, seq, n_seq, n_rep=1000):
    tri_idx = np.triu_indices(len(ec), 1)
    D = np.round(cdist(ec, ec), 2)[tri_idx]
    X = sorted(set(D))
    cov_mean = []
    for i in range(n_rep):
        cov = main_figs.get_seq_covariance(seq[np.random.choice(range(len(seq)), replace=False, size=n_seq)])[tri_idx]
        cov_mean.append([np.nanmean(cov[D==d]) for d in X])
    return np.array(cov_mean)


### Needs the results of the G-MeCh model (df2)
def plot_cov_vs_dist(df2, seq_all, ec):
    df = df2.loc[(df2.ID=='w2h1')&(df2.no=='293')&(df2.ang2==0.785398)]
    idx_list = [np.ones(len(df), bool), (df.gap>3.60)]
    fig, ax = plt.subplots(1,2,figsize=(8.0,3))
    fig.subplots_adjust(wspace=0.6)
    # Update the equilibrium coordinates with the location of
    # the chemical binding site
    ec = np.append(ec, ec[[1,6,2]], axis=0)
    tri_idx = np.triu_indices(len(ec), 1)
    D = np.round(cdist(ec, ec), 2)[tri_idx]

    lbls = ["All solutions", r"$\Delta \Delta G \geq 3.6$ kT"]
    col = [sns.color_palette("colorblind")[2], sns.color_palette("colorblind")[3]]
    X = sorted(set(D))[1:]

    xlim = [-0.05, 4.05]

    for i, idx in enumerate(idx_list):
        ax[i].plot(xlim, [0]*2, '-k', alpha=0.5)
        cov = main_figs.get_seq_covariance(np.array(list(df.loc[idx, 'pseq'])))[tri_idx]
        N = np.sum(idx)

        boot_path = PATH_DATA.joinpath(f"ran_cov_boot_{i}.npy")
        if boot_path.exists() and 1:
            boot_cov_ran = np.load(boot_path)
        else:
            boot_cov_ran = bootstrap_cov(ec, seq_all, N)
            np.save(boot_path, boot_cov_ran)

        Y1 = [np.nanmean(cov[D==d]) for d in X]
        Y1err = [np.nanstd(cov[D==d]) for d in X]
        
        Y2 = np.nanmean(boot_cov_ran, axis=0)[1:]
        Y2err = np.nanstd(boot_cov_ran, axis=0)[1:]

        ax[i].plot(X, Y1, '-o', color=col[i], label=f"{lbls[i]}\n$N_S={N}$", fillstyle='none')
#       ax[i].fill_between(X, Y2 - Y2err, Y2 + Y2err, color=col[i], alpha=0.2)
        ax[i].fill_between(X, Y2 - Y2err, Y2 + Y2err, color='grey', alpha=0.3)
        ax[i].set_xlim(*xlim)
        ax[i].legend(loc='lower right', frameon=False)
        ax[i].set_xlabel(r"$|r_i - r_j|$ / nm")
        ax[i].set_ylabel("Mean correlation")
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

    fig.savefig(PATH_FIG.joinpath(f"cov_vs_dist.pdf"), bbox_inches='tight')


#ef affinity_functional_forms(ds=5.5):
#   fig, ax = plt.subplots(3,2,figsize=(12,14))
#   ax = ax.reshape(ax.size)
#   fig.subplots_adjust(wspace=0.3, hspace=0.3)
#   a = 0.2
#   dt_arr = np.linspace(-1.0, 1.0, 1001)
#   alpha = np.zeros(dt_arr.size)
#   alpha[500] = 1
#   alpha[501:] = np.linspace(a*2, 0, 500)
#   alpha[250:500] = np.linspace(0, a, 250)

#   ax[0].plot(dt_arr, alpha)

#   Karr = 10.**np.linspace(1., 4., 1001)
#   ent = -ds * np.log10(Karr / 10000)
#   e = 20

#   ax[1].plot(np.log10(1/Karr), -ent)

#   def sig_fun(k, dt=0.01, e=8):
#       if isinstance(dt, float):
#           if dt == 0:
#               return np.ones(k.size)
#       f = (e / dt**2 / k)
#       return (f) / (f + 1)

#   idx = [350, 450, 500, 600, 700]
#   for i in idx:
#       sf = sig_fun(Karr, dt_arr[i], e)
#       lbl = r"$\Delta \theta_0$ = {0}".format(f"{dt_arr[i]:4.2f}")
#       ax[2].plot(np.log10(1/Karr), sf, label=lbl)
#       dg = e * (alpha[i] + (1-alpha[i]) * sf) - ent
#       ax[3].plot(np.log10(1/Karr), -np.max([np.zeros(dg.shape), dg], axis=0), label=lbl)

#   DT, K = np.meshgrid(dt_arr, Karr)
#   A, K = np.meshgrid(alpha, Karr)

#   DG = e * (alpha + (1 - alpha) * sig_fun(K, DT, e)).T - ent
#   im = ax[4].contourf(np.max([DG[:,::-1], np.zeros(DG.shape)], axis=0))
#   cbar = fig.colorbar(im, ax=ax[4])
#   cbar.ax.set_ylabel(r"Affinity, $\Delta G$ (kT)", rotation=270)
#   cbar.ax.get_yaxis().labelpad = 15
#   cbar.ax.set_yticklabels([-float(x.get_text()) if float(x.get_text())!=0 else 0.0 for x in cbar.ax.get_yticklabels()])

#   e = e / 2
#   DG = e * (alpha + (1 - alpha) * sig_fun(K, DT, e)).T - ent
#   im = ax[5].contourf(np.max([DG[:,::-1], np.zeros(DG.shape)], axis=0))
#   cbar = fig.colorbar(im, ax=ax[5])
#   cbar.ax.set_ylabel("Affinity, $\Delta G$ (kT)", rotation=270)
#   cbar.ax.get_yaxis().labelpad = 15
#   cbar.ax.set_yticklabels([-float(x.get_text()) if float(x.get_text())!=0 else 0.0 for x in cbar.ax.get_yticklabels()])

#   xlbls = [r"Shape mismatch, $\Delta \theta_0$"] + [r"Flexibility, $-\logK$"] * 5
#   ylbls = ["Zero-deformation binding affinity \n" + r"$\alpha~$ (kT)", r"$\Delta S$ (kT)",
#            r"$f_{deform}$", "Affinity, $\Delta G$ (kT)"] + [r"Shape mismatch, $\Delta \theta_0$"] * 5

#   for i, a in enumerate(ax):
#       a.set_xlabel(xlbls[i])
#       a.set_ylabel(ylbls[i])

#   xticks = [np.argmin(np.abs(Karr - 10.**x)) for x in range(1, 5)]
#   for i in [2,3]:
#       ax[i].legend(loc='best', frameon=False)
#   for i in [4,5]:
#       ax[i].set_xticks(xticks)
#       ax[i].set_xticklabels(range(-4,0,1))
#       ax[i].set_yticks([0, 500, 1000])
#       ax[i].set_yticklabels(range(-1, 2))

#   fs = 16
#   for i, b in enumerate('ABCDEF'):
#       ax[i].text( -0.20, 1.05, b, transform=ax[i].transAxes, fontsize=fs)
#   
#   fig.savefig(PATH_FIG.joinpath(f"func_form.pdf"), bbox_inches='tight')


### Needs the results of the MeCh
def spec_aff_def(base=PATH_BASE, ds=5.5, e=8, it=12, fig='', ax='', ft=14, km=226, ks=10000):
    if isinstance(ax, str):
        fig = plt.figure(figsize=(12, 5))
        gs = GridSpec(3,3, height_ratios=[.2, 1, 1])
        ax = [fig.add_subplot(gs[1:,0])] + \
             [fig.add_subplot(gs[i,1]) for i in [1,2]] + \
             [fig.add_subplot(gs[1:,2])]
        cax = [fig.add_subplot(gs[0,i]) for i in [0,1,2]]
        fig.subplots_adjust(hspace=0.05, wspace=0.4)

#   cmap = Oslo_20.mpl_colormap.reversed()
    cmap1 = YlOrBr_9.mpl_colormap
    cmap2 = YlGn_9.mpl_colormap    
    cmap3 = Purples_9.mpl_colormap    

    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]

    Earr= np.arange(0.25, 25.25, 0.25)
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    i = np.argmin(np.abs(Earr - e))
    ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"{i:03d}", "energy.out")) for i in range(km)])
    ener = utils.update_energy_def(ener[:km,:], ang, Earr[i])#/ e
    edef = ener[:,:,1]
    ener = ener[:,:,0]

    gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
    gap2 = np.abs(ener[:,14:24] - ener[:,12:22])
    ener = ener[:,12:]
    edef = edef[:,12:]
#   ang = ang[12:]

    ent = np.array([-ds * np.log10(Kmax/ks)] * ener.shape[1]).T
    fener = np.abs(np.min([ener + ent, np.zeros(ener.shape)], axis=0))
    im = []
    im.append(ax[0].contourf(fener.T[:,::-1], cmap=cmap1))

    im.append(ax[1].contourf(gap1.T[:,::-1], cmap=cmap2, vmin=0, vmax=e))
    im.append(ax[2].contourf(gap2.T[:,::-1], cmap=cmap2, vmin=0, vmax=e))

    im.append(ax[3].contourf(edef.T[:,::-1], cmap=cmap3))

    # Plot line of optimal mismatch for given flexbility
    ce1 = 'k'
    ce2 = 'k'
    X = np.arange(km)
    Y = np.argmax(gap1, axis=1)[::-1]
    main_figs.plot_optimal_line(ax[1], X, Y, 9, ce1)
    [main_figs.plot_optimal_line(ax[i], X, Y, 9, ce1, offset=11) for i in [0,3]]
    [main_figs.plot_optimal_line(ax[i], X, Y, 9, ce2, offset=13, p=':') for i in [0,3]]

    Y = np.argmax(gap2, axis=1)[::-1]
    main_figs.plot_optimal_line(ax[2], X, Y, 0, ce1)
    [main_figs.plot_optimal_line(ax[i], X, Y, 0, ce1, offset=2) for i in [0,3]]
    [main_figs.plot_optimal_line(ax[i], X, Y, 0, ce2, offset=0, p=':') for i in [0,3]]

    # Plot line of optimal mismatch for given flexbility for higher energy (e=16)
    i = 79
    ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
    ener = main_figs.update_energy(ener[:km,:], ang, Earr[i])#/ e
    gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
    gap2 = np.abs(ener[:,14:24] - ener[:,12:22])

#   ce2 = 'purple'
#   Y = np.argmax(gap1, axis=1)[::-1]
#   plot_optimal_line(ax[1], X, Y, 9, ce2)

#   Y = np.argmax(gap2, axis=1)[::-1]
#   plot_optimal_line(ax[2], X, Y, 0, ce2)

    # Add colorbar
    clbl = [r"$\Delta G_C$ (kT)", r"$\Delta \Delta G$ (kT)", r"$E_{def}$ (kT)"]
    cticks = [np.arange(0, 30, 8), np.arange(0, 8, 2), np.arange(0, 8, 2)]
    for i, j in enumerate([0,2,3]):
        cbar = fig.colorbar(im[j], cax[i], pad=0.2, orientation='horizontal', ticks=cticks[i])
        cbar.set_label(clbl[i], rotation=0, fontsize=ft-4)
        cbar.ax.get_yaxis().labelpad = 20
        cax[i].xaxis.set_ticks_position('top')
        cax[i].xaxis.set_label_position('top')
    cax[0].set_xticklabels([-float(x.get_text()) if float(x.get_text())!=0 else 0.0 for x in cax[0].get_xticklabels()])
        
    # Add labels, ticks, etc.
    ylbl = [r"Shape mismatch, $\Delta \theta_{0} ~(^{{\circ}})$"] + \
           [r"$\Delta \theta_{0} ~(^{{\circ}})$"] * 3
    xticks = [np.argmin(np.abs(Kmax - 10.**x)) for x in range(1, 5)]
    ylim = [22, 9, 9, 22]
    for i in range(4):
        ax[i].set_xlabel(r"Flexibility, $-\log_{10} K$")
        ax[i].set_ylabel(ylbl[i])
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(range(-4,0,1))
        ax[i].set_ylim(0, ylim[i])
    ax[1].set_xticks([])
    ax[1].set_xlabel("")
    ax[0].set_yticks(np.arange(3, 22, 4))
    ax[0].set_yticklabels(np.arange(-20, 30, 10))
    ax[1].set_yticks(np.arange(0, 12, 4))
    ax[1].set_yticklabels(np.arange(0, 30, 10))
    ax[1].set_yticklabels(["", "10", "20"])
    ax[2].set_yticks(np.arange(1, 12, 4))
    ax[2].set_yticklabels(["-20", "-10", ""])
    ax[3].set_yticks(np.arange(3, 22, 4))
    ax[3].set_yticklabels(np.arange(-20, 30, 10))

    fig.savefig(PATH_FIG.joinpath(f"spec_aff_def.pdf"), bbox_inches='tight')


### Needs the results of the MeCh
def compare_gradients(e=8, ds=5.5, km=226):
    ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter12", f"{i:03d}", "energy.out")) for i in range(km)])
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter12", f"Template", "inputs", "ligands.dat"))[:,3]
    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ener = utils.update_energy_def(ener, ang, e)
    ent = np.array([ds*np.log10(10000/Kmax)] * 35).T
    fener = ener[:,:,0] + ent
    f1 = fener[:,14:24]
    f2 = fener[:,23:-2]
    gap1 = fener[:,14:24] - fener[:,12:22]
    gap2 = fener[:,23:-2] - fener[:,25:]
    def1 = ener[:,14:24,1]
    def2 = ener[:,23:-2,1]
#   def1 = ener[:,12:22,1]
#   def2 = ener[:,25:,1]

#   return fener, f1, f2, gap1, gap2, def1, def2

#   fig, ax = plt.subplots(2,4)
    fig, ax = plt.subplots(5,2)
    ax[0,0].contourf(f1)
    ax[0,1].contourf(f2)
    ax[1,0].contourf(gap1)
    ax[1,1].contourf(gap2)
    ax[2,0].contourf(def1)
    ax[2,1].contourf(def2)

    for a, f, g, d in zip(ax[3:].T, [f1, f2], [gap1, gap2], [def1, def2]):
#   for f, g, d in zip([f1, f2], [gap1, gap2], [def1, def2]):
        fx, fy = np.gradient(f)
        gx, gy = np.gradient(g)
        dx, dy = np.gradient(d)

        gcut = 10
        dfcut = 0.1
        dgcut = 0.1

#       idx = ((fx < 0) & (gx < 0) & (dx > 0)) | ((fx > 0) & (gx > 0) & (dx < 0))
#       idx = ((fx < 0) & (gx < 0) & (dx < 0)) | ((fx > 0) & (gx > 0) & (dx > 0))
        idx = ((fx < 0) & (gx < 0) & (dx < 0)) | ((fx > 0) & (gx > 0) & (dx > 0)) | ((fy < 0) & (gy < 0) & (dy < 0)) | ((fy > 0) & (gy > 0) & (dy > 0))

#       idx = ((fx < 0) & (gx < 0)) | ((fx > 0) & (gx > 0)) | ((fy < 0) & (gy < 0)) | ((fy > 0) & (gy > 0))
#       idx = ((fx < 0) & (gx < 0)) | ((fx > 0) & (gx > 0))
#       idx = ((fy < 0) & (gy < 0)) | ((fy > 0) & (gy > 0))
#       idx = ((fx < 0) & (gx < 0)) | ((fx > 0) & (gx > 0))
        im = a[0].imshow(idx.astype(int), aspect=5/226, cmap='gray')
        cbar = fig.colorbar(im, ax=a[0])

#       idx = ((fy < 0) & (gy < 0) & (dy > 0)) | ((fy > 0) & (gy > 0) & (dy < 0))
#       idx = ((fy < 0) & (gy < 0) & (dy < 0)) | ((fy > 0) & (gy > 0) & (dy > 0))
#       idx = ((fx < 0) & (gx < 0) & (dx < 0)) | ((fx > 0) & (gx > 0) & (dx > 0)) | ((fy < 0) & (gy < 0) & (dy < 0)) | ((fy > 0) & (gy > 0) & (dy > 0))

        idx = ((fx < 0) & (gx < 0) & (dx > 0)) | ((fx > 0) & (gx > 0) & (dx < 0)) | ((fy < 0) & (gy < 0) & (dy > 0)) | ((fy > 0) & (gy > 0) & (dy < 0))
#       idx = ((fx < 0) & (gx < 0) & (dx > 0)) | ((fx > 0) & (gx > 0) & (dx < 0))
#       idx = ((fy < 0) & (gy < 0) & (dy > 0)) | ((fy > 0) & (gy > 0) & (dy < 0))
#       idx = ((fy < 0) & (gy < 0)) | ((fy > 0) & (gy > 0))
        im = a[1].imshow(idx.astype(int), aspect=5/226, cmap='gray')
        cbar = fig.colorbar(im, ax=a[1])
#       sns.distplot(dx[idx] * -np.sign(fx[idx]))

        a[0].invert_yaxis()
        a[1].invert_yaxis()


### Needs the results of the MeCh
def spec_aff_de(base=PATH_BASE, de=1, ds=5.5, e=8, it=12, fig='', ax='', ft=14, km=226, ks=10000):
    if isinstance(ax, str):
        fig = plt.figure(figsize=(7, 4.8))
        gs = GridSpec(2,3, height_ratios=[.2, 1])
        ax = [fig.add_subplot(gs[1,i]) for i in range(3)]
        cax = [fig.add_subplot(gs[0,i]) for i in range(3)]
        fig.subplots_adjust(hspace=0.05, wspace=0.4)

#   cmap = Oslo_20.mpl_colormap.reversed()
    cmap1 = YlOrBr_9.mpl_colormap
    cmap2 = YlGn_9.mpl_colormap    

    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]

    Earr= np.arange(0.25, 25.25, 0.25)
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    i = np.argmin(np.abs(Earr - e))
    j = np.argmin(np.abs(Earr - (e-de)))

    ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
    ener1 = main_figs.update_energy(ener[:km,:], ang, Earr[i])#/ e

    ener = np.load(base.joinpath(f"{j:03d}", "energy.dat.npy"))[:,:,0]
    ener2 = main_figs.update_energy(ener[:km,:], ang, Earr[j])#/ e

    gap = ener2[:,14:-2] - ener1[:,14:-2]
    ener1 = ener1[:,12:]
    ener2 = ener2[:,12:]
#   ang = ang[12:]

    ent = np.array([1.5 - 2 * np.log10(Kmax/ks)] * ener1.shape[1]).T + 15
#   ent = np.ones(ent.shape) * -5.5
    fener1 = np.abs(np.min([ener1 + ent, np.zeros(ener1.shape)], axis=0))
    fener2 = np.abs(np.min([ener2 + ent, np.zeros(ener2.shape)], axis=0))
    im = []
    im.append(ax[0].contourf(fener1.T[:,::-1], cmap=cmap1))
    im.append(ax[1].contourf(fener2.T[:,::-1], cmap=cmap1))

    im.append(ax[2].contourf(gap.T[:,::-1], cmap=cmap2, vmin=0, vmax=de*3))

    # Plot line of optimal mismatch for given flexbility
#   ce1 = 'k'
#   X = np.arange(km)
#   Y = np.argmax(gap, axis=1)[::-1]
#   main_figs.plot_optimal_line(ax[0], X, Y, 9, ce1, offset=11)
#   main_figs.plot_optimal_line(ax[1], X, Y, 9, ce1)

#   Y = np.argmax(gap2, axis=1)[::-1]
#   main_figs.plot_optimal_line(ax[0], X, Y, 0, ce1, offset=2)
#   main_figs.plot_optimal_line(ax[2], X, Y, 0, ce1)


    # Add colorbar
    clbl = [r"$\Delta G_C$ (kT)", r"$\Delta G_C$ (kT)", r"$\Delta \Delta G$ (kT)"]
    cticks = [np.arange(0, 30, 8), np.arange(0, 30, 8), np.arange(0, 4*de, de)]
    for i in range(3):
        cbar = fig.colorbar(im[i], cax[i], pad=0.2, orientation='horizontal', ticks=cticks[i])
        cbar.set_label(clbl[i], rotation=0, fontsize=ft-4)
        cbar.ax.get_yaxis().labelpad = 20
        cax[i].xaxis.set_ticks_position('top')
        cax[i].xaxis.set_label_position('top')
    cax[0].set_xticklabels([str(-1 * int(x.get_text())) for x in cax[0].get_xticklabels()])
        
    # Add labels, ticks, etc.
    ylbl = [r"Shape mismatch, $\Delta \theta_{0} ~(^{{\circ}})$",
            r"$\Delta \theta_{0} ~(^{{\circ}})$", r"$\Delta \theta_{0} ~(^{{\circ}})$"]
    xticks = [np.argmin(np.abs(Kmax - 10.**x)) for x in range(1, 5)]
    ylim = [22, 22, 18]
    for i in range(3):
        ax[i].set_xlabel(r"Flexibility, $-\log_{10} K$")
        ax[i].set_ylabel(ylbl[i])
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(range(-4,0,1))
        ax[i].set_ylim(0, ylim[i])
    ax[1].set_xticks([])
    ax[1].set_xlabel("")
    ax[0].set_yticks(np.arange(3, 22, 4))
    ax[0].set_yticklabels(np.arange(-20, 30, 10))
    ax[1].set_yticks(np.arange(0, 12, 4))
    ax[1].set_yticklabels(np.arange(0, 30, 10))
    ax[1].set_yticklabels(["", "10", "20"])
    ax[2].set_yticks(np.arange(1, 12, 4))
    ax[2].set_yticklabels(["-20", "-10", ""])



### Needs the sequences used in the G-MeCh model
def entropy_scaling(seq2):
    fig, ax = plt.subplots()

    params = {'text.usetex': True,
              "text.latex.preamble": [r"\usepackage{amsmath}"]}
    


    c0 = utils.load_equil_coord()
    adj = rec_utils.get_adj_from_coord(c0)
    idx_break =[(1,2), (1,6), (2,6)]
    adj = rec_utils.break_links(adj, idx_break)
    Krat = 10.**np.linspace(-4, 4, 101)

    ent = np.array([GE.entropy_diff(c0, adj, np.array([[1,1],[1,1]]), seq2[0], idx0=[1,6,2], kstiff=k) for k in Krat])
    ax.plot(Krat, ent, 'o', fillstyle='none')

    idx = Krat > 10
    slope, intercept = linregress(np.log(Krat[idx]), ent[idx])[:2]
    X1 = 10.**np.linspace(-5, 5)
    Y1 = slope * np.log(X1) + intercept
    ax.plot(X1, Y1, '--k')
    print(slope, intercept)

    fs = 12
    
    s = r"{0} + {1} $\ln{{K_{{\Lambda}}/K}}$".format(round(intercept, 1), round(slope, 1))
    ax.text(0.54, 0.61, s, transform=ax.transAxes, fontsize=fs, rotation=40)

    idx = Krat < 0.1
    slope, intercept = linregress(np.log(Krat[idx]), ent[idx])[:2]
    Y2 = slope * np.log(X1) + intercept
    ax.plot(X1, Y2, ':k')
    print(slope, intercept)

    s = r"{0} + {1} $\ln{{K_{{\Lambda}}/K}}$".format(round(intercept, 1), round(slope, 1))
    ax.text(0.04, 0.20, s, transform=ax.transAxes, fontsize=fs, rotation=30)

    ax.set_xscale('log')
    ax.set_xlabel(r"$K_{{\Lambda}}$ / $K$")
    ax.set_ylabel(r"$-T\Delta S_{conf}$ / kT")

    fig.savefig(PATH_FIG.joinpath(f"si8.pdf"), bbox_inches='tight')


### Needs the results of evolve_function.py
def fitness_landscape_minima():
    fig, ax = plt.subplots(5,3,figsize=(10,8))
    fig.subplots_adjust(hspace=0.5,wspace=0.4)
    base = Path("/data/Jmcbride/RecProt/Data/EvolveSim")
    path_list = [base.joinpath(f"{f}.pkl") for f in ["OpenOpen_minima", "LocalStruct12_minima", "LocalStruct12_minima_2"]]
    res = [pickle.load(open(f, 'rb')) for f in path_list]

    fc_list1 = np.round(np.arange(0.96, 0.78, -0.04), 2)
    fc_list2 = np.round(np.arange(0.4,-0.1, -0.1), 1)

    X = np.arange(20, 110, 10)
    for i in range(5):
        ax[i,0].plot(X, res[0]['results'][fc_list1[i]][3][2:-2][::-1], '-s', fillstyle='none', label=r"$f_{cut} = $" + f"{fc_list1[i]}", c=sns.color_palette()[0])
        ax[i,1].plot(X, res[1]['results'][fc_list1[i]][3][2:-2][::-1], '-s', fillstyle='none', label=r"$f_{cut} = $" + f"{fc_list1[i]}", c=sns.color_palette()[1])
        ax[i,2].plot(X, res[2]['results'][fc_list2[i]][3][2:-2][::-1], '-s', fillstyle='none', label=r"$f_{cut} = $" + f"{fc_list2[i]}", c=sns.color_palette()[2])

        for j in range(3):
            ax[i,j].legend(loc='upper right', frameon=False)
            ax[i,j].plot([X[0], X[-1]], [1,1], ':k')
            ax[i,j].spines['top'].set_visible(False)
            ax[i,j].spines['right'].set_visible(False)
            ax[i,j].set_ylim(-1, 27)
            ax[i,j].set_xticks(np.arange(20, 120, 20))
    ax[2,0].set_ylabel("Number of local minima in fitness landscape")

    ttls = ["G-MeCh; one non-cognate", "G-MeCh-S; one non-cognate", "G-MeCh-S; two non-cognates"]

    for j in range(3):
        ax[0,j].set_title(ttls[j], loc='left', fontsize=11)
        ax[-1,j].set_xlabel(r"$\theta_C~(^{{\circ}})$")
    ax[2,0].set_ylabel("Number of basins of attraction in fitness landscape")

    fig.savefig(PATH_FIG.joinpath(f"si9.pdf"), bbox_inches='tight')


def get_theta(c, idx):
    v1 = c[idx[1]] - c[idx[0]]
    v2 = c[idx[1]] - c[idx[2]]
    return np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2))


### Needs the results of the G-MeCh model
def plot_keff_vs_gap(Keff='', G='', Eex=''):
    emat = np.array([3.875, 6.0])
    if isinstance(Keff, str) or isinstance(G, str) or isinstance(Eex, str):
        seq = np.array([list(x) for x in np.load('../Data/Sequences/seq_str_2.npy')]).astype(int)
        seq_ener = emat[seq[:,13:]-1]
        Eex = np.min(seq_ener[:,[0,2]], axis=1)

        C = np.loadtxt('/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/w2h1/273/config.out').reshape(65536,13,13,2)
        E = np.loadtxt('/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/w2h1/273/energy.out').reshape(2**16,13,4)
        F = np.load('/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/w2h1/273/free_energy.npy').reshape(2**16,13)
        G = np.abs(F[:,:-1] - F[:,1:])

        T = np.array([[get_theta(C[i,j], [1,6,2]) for j in range(13)] for i in range(65536)])
        R0 = np.sin(T - np.pi / 3)
        Keff = E[:,:,1] / R0**2
    

    ang = (90 - np.arange(30, 95, 5)) * 2


    fig, ax = plt.subplots(2,2,figsize=(10,10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    col = Berlin_11.mpl_colors
    col = sns.color_palette()
    col = Acton_6.mpl_colors
    ms = 6
    al = 0.9
    ttl2 = ['Open pocket', 'Narrow pocket']

    for k, istart in enumerate([0, 7]):
        idx = np.arange(istart, len(Keff), 8)[::20]
        for j, i in enumerate(range(5)[::-1]):
            lbl = r"$\Delta \theta_0 = {0}$".format(int(round(60 - ang[i+1])))
            ax[k,1].plot(-np.log10(Keff[idx,i+1]), G[idx,i]/Eex[idx], 'o', label=lbl, c=col[j], ms=ms, alpha=al, mec='k', mew=0.2)

        for j, i in enumerate(range(7,12)):
            lbl = r"$\Delta \theta_0 = {0}$".format(int(round(60 - ang[i])))
            ax[k,0].plot(-np.log10(Keff[idx,i]), G[idx,i]/Eex[idx], 'o', label=lbl, c=col[j], ms=ms, alpha=al, mec='k', mew=0.2)

        ax[k,0].set_title(r"$\epsilon = {0}$ kT, {1}".format(emat[k], ttl2[0]), loc='center')
        ax[k,1].set_title(r"$\epsilon = {0}$ kT, {1}".format(emat[k], ttl2[1]), loc='center')

    fs = 14
    for i, a in enumerate(ax.ravel()):
        a.legend(loc='best', frameon=False)
        a.set_xlim(-3, -0.5)
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_xlabel(r"Flexibility, $-\log_{10}~K_{eff}$ (kT/nm$^2$)")
        a.set_ylabel(r"$\Delta\Delta G $ (kT)")
        a.text( -0.10, 1.03, 'ABCD'[i], transform=a.transAxes, fontsize=fs)

    fig.savefig(PATH_FIG.joinpath(f"si10.pdf"), bbox_inches='tight')


### Needs the results of the MeCh model
def spec_vs_e(base=PATH_BASE, ds=0.0, e=8, it=12, fig='', ax='', cax='', ft=14, km=226, ks=10000):
    if isinstance(ax, str):
        fig = plt.figure(figsize=(8, 5))
        gs = GridSpec(3,3, height_ratios=[.2, 1, 1])
        ax = [[fig.add_subplot(gs[i,j]) for i in [1,2]] for j in [0,1,2]]
        cax = [fig.add_subplot(gs[0,i]) for i in [0,1,2]]
        fig.subplots_adjust(hspace=0.05, wspace=0.5)

    cmap = YlGn_9.mpl_colormap    
    ncont = 16

    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]

    Earr= np.arange(0.25, 25.25, 0.25)
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    e_print = [2, 6, 12]
    im = []


    def plot_optimal_line(ax, X, Y, cut, col, offset=0, p='--', d=2):
        idx = range(np.where(Y==cut)[0][0] + 10)
        X, Y = X[idx], Y[idx]
#       Yp = np.poly1d(np.polyfit(X, Y, d))(X)
#       ax.plot(X, Yp + offset, p, color=col)
        ax.plot(X, Y, p, color=col)
        

    for a, e in zip(ax, e_print):

        i = np.argmin(np.abs(Earr - e))
        ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
        ener = main_figs.update_energy(ener[:km,:], ang, Earr[i])#/ e

        gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
        gap2 = np.abs(ener[:,14:24] - ener[:,12:22])
        ener = ener[:,12:]

        ent = np.array([1.5 - 2 * np.log(Kmax/ks)] * ener.shape[1]).T + ds
        fener = np.abs(np.min([ener + ent, np.zeros(ener.shape)], axis=0))
        im.append(a[0].contourf(gap1.T[:,::-1], cmap=cmap, vmin=0, vmax=e, levels=ncont))
        im.append(a[1].contourf(gap2.T[:,::-1], cmap=cmap, vmin=0, vmax=e, levels=ncont))

        ce1 = 'k'
        X = np.arange(km)
        Y = np.argmax(gap1, axis=1)[::-1]
        plot_optimal_line(a[0], X, Y, 9, ce1, d=2)

        Y = np.argmax(gap2, axis=1)[::-1]
        plot_optimal_line(a[1], X, Y, 0, ce1, d=2)

    print(len(im))

    # Add colorbar
    cticks = [np.arange(0, 3, 1), np.arange(0, 8, 2), np.arange(0, 14, 4)]
    for i, j in enumerate([1,3,5]):
        cbar = fig.colorbar(im[j], cax[i], pad=0.2, orientation='horizontal', ticks=cticks[i])
        cbar.set_label(r"$\Delta \Delta G$ (kT)", rotation=0, fontsize=ft-4)
        cbar.ax.get_yaxis().labelpad = 20
        cax[i].xaxis.set_ticks_position('top')
        cax[i].xaxis.set_label_position('top')
        
    # Add labels, ticks, etc.
    fs = 12
    xticks = [np.argmin(np.abs(Kmax - 10.**x)) for x in range(1, 5)]
    for i in range(3):
        ax[i][0].set_yticklabels(np.arange(0, 50, 20))
        ax[i][0].set_yticklabels(["0", "20", "40"])
        ax[i][1].set_yticks(np.arange(1, 12, 4))
        ax[i][1].set_yticklabels(["-40", "-20", ""])

        ax[i][0].set_ylabel(r"$\Delta \theta_{0} ~(^{{\circ}})$")
        ax[i][1].set_ylabel(r"$\Delta \theta_{0} ~(^{{\circ}})$")
        ax[i][0].set_xlabel("")
        ax[i][1].set_xlabel(r"Flexibility, $-\log_{10} K$")
        ax[i][0].set_xticks([])
        ax[i][1].set_xticks(xticks)
        ax[i][1].set_xticklabels(range(-4,0,1))
        cax[i].text( -0.20, 1.05, 'ABCD'[i], transform=cax[i].transAxes, fontsize=fs)


    fig.savefig(PATH_FIG.joinpath(f"si11.pdf"), bbox_inches='tight')






