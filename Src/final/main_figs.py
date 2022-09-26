"""
Code for the figures in the main manuscript.
Will not run without the correct data, or without correct paths to data.

Fig 1, 2 & 5 require the main FORTRAN code to be run on the no-sequence MeCh model;
    in practice, this is done by using a sequence with only one amino acid type.

Fig 3 requires the G-MeCh model to be run

Fig 4 requires the output of "evolve_function.py" using as input
    the results of the G-MeCh-S model for many protein shapes

Fig 6 does not need any data

"""
import os
from pathlib import Path
import pickle

import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from palettable.colorbrewer.qualitative import Paired_12
from palettable.colorbrewer.sequential import Purples_9, YlGn_9, YlOrBr_9
from palettable.scientific.sequential import Acton_20, Oslo_20, Nuuk_20, LaJolla_20, LaPaz_20
from palettable.scientific.diverging import Berlin_20, Roma_20, Berlin_10, Vik_20, Vik_6
import pandas as pd
from PIL import Image
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
from scipy.stats import linregress, pearsonr
import seaborn as sns

import local_utils

PATH_FIG = Path("/data/Jmcbride/RecProt/Figures/Paper")
PATH_DATA = Path("/data/Jmcbride/RecProt/Figures/Paper/Data")
PATH_BASE = Path("/molahome/jmcbride/RecProt")



########################################################
### FIG 1
###

def fig1():
    fig = plt.figure(figsize=(14,5))
    gs = GridSpec(2,2, width_ratios=[1, 1.2])
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    ax = [fig.add_subplot(gs[i,j]) for j in range(2) for i in range(2)]

    ax[0].set_axis_off()

    optimal_mismatch_flex(fig=fig, ax=ax[1:], ft=11)

    for a in ax:
        for direction in ['right', 'top']:
            a.spines[direction].set_visible(False)

    fs = 16
    for i, b in enumerate('ABCD'):
        ax[i].text( -0.20, 1.05, b, transform=ax[i].transAxes, fontsize=fs)

    fig.savefig(PATH_FIG.joinpath("fig1.pdf"), bbox_inches='tight')



def update_energy(ener, ang, e):
    theta = np.pi/3
    idx = ang < theta
    ener[:,idx] = np.min([ener[:,idx], np.zeros((ener.shape[0], idx.sum())) - e], axis=0)
    idx = ang == theta
    ener[:,idx] = np.min([ener[:,idx], np.zeros((ener.shape[0], idx.sum())) - e*3], axis=0)
    idx = ang > theta
    ener[:,idx] = np.min([ener[:,idx], np.zeros((ener.shape[0], idx.sum())) - e*2], axis=0)
    return ener


def get_rigid_binding_energy(angle, e, a0=np.pi/3, sigma=0.3):
    chord = 2 * np.sin((a0 - angle)/2)
    if angle < a0:
        return np.sum([- e * np.exp(-(r * chord)**2/sigma**2) for r in [0, 0, 1]])
    elif abs(angle - a0) < 0.001:
        return -e * 3
    elif angle > a0:
        return -e


def optimal_mismatch_vs_mismatch(fig='', ax='', sigma=0.3):
    if isinstance(ax, str):
        fig, ax = plt.subplots(2, 1, figsize=(10,8))

    dt = np.pi/180 * 1
    thetaB = np.arange(np.pi/3, dt, -np.pi/180)
    dtheta = 180/np.pi * (thetaB[:-1][::-1] - thetaB[-1])

    dE = np.exp(-(2*np.sin(np.pi/3 - thetaB)/2)**2 / sigma**2)
    opt_dt = 60 - 180/np.pi * np.array([thetaB[np.argmax(dE[:-i] - dE[i:])] for i in range(1, thetaB.size)])
    max_gap = np.array([np.max(dE[:-i] - dE[i:]) for i in range(1, thetaB.size)])

    ax[0].plot(dtheta, opt_dt)
    ax[1].plot(dtheta, max_gap)

    ax[0].set_xlabel(r"$\theta_L - \theta_S$")
    ax[1].set_xlabel(r"$\theta_L - \theta_S$")
    ax[0].set_ylabel(r"Optimal $\theta_0 - \theta_L$")
    ax[1].set_ylabel(r"Max $\Delta \Delta G$")


def optimal_mismatch_flex(base=PATH_BASE, e=8, it=12, fig='', ax='', ft=14, km=225):
    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]
    ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"{i:03d}", "energy.out")) for i in range(km)])[:,:,0]
#   ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"AllE", f"023", f"{i:03d}", "energy.out")) for i in range(300)])[:,:,0]
    ener = update_energy(ener, ang, e) / e
    gap = ener[:,:-2] - ener[:,2:]

    if isinstance(ax, str):
        fig, ax = plt.subplots(3, 1, figsize=(10,10))
        fig.subplots_adjust(hspace=0.3)

    ang2 = np.arange(0, np.pi/2, np.pi/180)
#   X2 = 2 * (90 - ang2*180/np.pi)
    X2 = ang2 * 180 / np.pi
    E2 = np.array([get_rigid_binding_energy(a, 1) for a in ang2])
    for i, dtheta in enumerate([5, 10, 20]):
        x, y = X2[dtheta:], E2[dtheta:] - E2[:-dtheta]
        idx = (x >=30) & (x <= 90)
        c = sns.color_palette()[i]
        ax[0].plot(x[idx], y[idx], label=r"$\Delta \theta_{{LS}} = {0}^{{\circ}}$".format(dtheta), c=c)
        imin = np.argmin(y[idx])
        ax[0].plot([x[idx][imin]]*2, [-1.2, y[idx][imin]], ':', c=c)
    ax[0].set_yticks(np.arange(-1, 2.5, 1))
    ax[0].set_xlabel(r"Larger ligand, $\theta_L$ (degrees)")
#   ax[0].set_ylabel("Binding Energy gap\n $\\Delta G_B - \\Delta G_A$", fontsize=ft)
    ax[0].set_ylabel(r"$\Delta \Delta G_{LS}$ ($\epsilon$ kT)", fontsize=ft)
    ax[0].legend(loc='upper left', frameon=False)
    ax[0].set_ylim(-1.2, 2.2)

    X = 2 * (90 - ang[:-2]*180/np.pi)
    theta = np.pi/3
    j0 = np.argmin(np.abs(ang - theta))
    idx_list1 = sorted(set(list(range(1,gap.shape[1],4)) + [j0 - 2, j0]))
    for i, j in enumerate(idx_list1):
        if j == j0 - 2:
            ax[1].plot(Kmax, gap[:,j], '--k', lw=2)
        elif j == j0:
            ax[1].plot(Kmax, gap[:,j], ':k', lw=2)
        else:
            ax[1].plot(Kmax, gap[:,j], color=Berlin_20.mpl_colormap(j/(gap.shape[1]-1)))
#           ax[1].plot(Kmax, gap[:,j], color=Vik_20.mpl_colormap(j/(gap.shape[1]-1)))

    idx_list2 = range(0, gap.shape[0], 20)
    for i, j in enumerate(idx_list2):
        if i == 0:
            ax[2].plot(X+5, gap[j], ':k', lw=2)
        else:
#           k = j/(gap.shape[0]-1)
#           k = k / (1-k)
            k = np.log(Kmax[j]) / np.log(Kmax[-1])
            ax[2].plot(X, gap[j], color=Vik_20.mpl_colormap(k))
    ax[2].plot(X, gap[-1], '-k', lw=2)

    norm = Normalize(vmin=X[idx_list1[-1]], vmax=X[idx_list1[0]])
    cax1 = fig.colorbar(ScalarMappable(norm=norm, cmap=Berlin_20.mpl_colormap), ax=ax[1])
    cax1.ax.get_yaxis().labelpad = 15
    cax1.ax.set_ylabel(r"$\theta_L$", rotation=0, fontsize=ft-1)


    norm = LogNorm(vmin=Kmax[idx_list2[0]], vmax=Kmax[idx_list2[-1]])
    cax2 = fig.colorbar(ScalarMappable(norm=norm, cmap=Vik_20.mpl_colormap), ax=ax[2])
    cax2.ax.get_yaxis().labelpad = 15
    cax2.ax.set_ylabel(r"$K$ (kT/nm$^2$)", rotation=270, fontsize=ft-1)


    ax[1].set_ylim(-0.505, 0.75)
    ax[1].set_xlim(7, 50000)
    ax[1].set_xscale('log')
    ax[1].set_xlabel(r"Spring Constant, $K$ (kT/nm$^{2}$)", fontsize=ft)
#   ax[1].set_ylabel("Binding Energy gap\n $\\Delta G_B - \\Delta G_A$", fontsize=ft)
    ax[1].set_ylabel(r"$\Delta \Delta G_{LS}$ ($\epsilon$ kT)", fontsize=ft)
    ax[2].set_xlabel(r"$\theta_L$ (degrees)", fontsize=ft)
#   ax[2].set_ylabel("Binding Energy gap\n $\\Delta G_B - \\Delta G_A$", fontsize=ft)
    ax[2].set_ylabel(r"$\Delta \Delta G_{LS}$ ($\epsilon$ kT)", fontsize=ft)

    ax[2].set_xticks(np.arange(30, 180, 30))

    for a in ax:
        for tick in a.xaxis.get_major_ticks():
            tick.label.set_fontsize(ft)
        for tick in a.yaxis.get_major_ticks():
            tick.label.set_fontsize(ft)

    ax[1].annotate(r"$\theta_L$", xy=(20,-0.4), xycoords='data', xytext=(20,  0.5), horizontalalignment='left',
                   arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=ft)

    i0 = np.argmin(np.abs(Kmax - 1000))
    ax[1].annotate(r"$\theta_S=60$", xy=(Kmax[i0+25], gap[i0+25,j0-2]), xycoords='data', xytext=(45000, 0.2), horizontalalignment='right', verticalalignment='bottom',
                   arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=ft)

    i0 = np.argmin(np.abs(Kmax - 10000))
    ax[1].annotate(r"$\theta_L=60$", xy=(Kmax[i0-25], gap[i0,j0]), xycoords='data', xytext=(45000, -0.25), horizontalalignment='right', verticalalignment='top',
                   arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=ft)

    ax[2].annotate(r"$K$", xy=(60,0.5), xycoords='data', xytext=(80, -0.4), horizontalalignment='left',
                   arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=ft)

#   fig.savefig(PATH_FIG.joinpath("opt_mismatch_flex.pdf"), bbox_inches='tight')


########################################################
### FIG 2
###


def fig2(i0=3):
    fig = plt.figure(figsize=(15,7))
    gs = GridSpec(22,6, width_ratios=[1,0.3,1.5,0.25,1,0.05])
    fig.subplots_adjust(hspace=2.0, wspace=0)
    ax = [fig.add_subplot(gs[1:10,0]), fig.add_subplot(gs[13:,0])] + \
         [fig.add_subplot(gs[2:12,2])] + \
         [fig.add_subplot(gs[i*6+i*2:(i+1)*6+i*2,4]) for i in [0,1,2]] + \
         [fig.add_subplot(gs[15:,2])]

    cax = [fig.add_subplot(gs[:2,2]), fig.add_subplot(gs[2:-2,5])]

    ft = 12
    for a in ax:
        for tick in a.xaxis.get_major_ticks():
            tick.label.set_fontsize(ft)
        for tick in a.yaxis.get_major_ticks():
            tick.label.set_fontsize(ft)


    entropy_effect(fig=fig, ax=ax[:2], ft=ft)
    plot_rec_space_noseq(fig=fig, ax=ax[2:], ft=ft, cax=cax, i0=i0)
    specificity_vs_affinity(fig=fig, ax=ax[-1], ft=ft, i0=i0)

    fs = 16
    ax[0].text(-0.18, 1.30, 'A', transform=ax[0].transAxes, fontsize=fs)
    cax[0].text(-0.15, 1.35, 'B', transform=cax[0].transAxes, fontsize=fs)
    ax[3].text( -0.25, 1.10, 'C', transform=ax[3].transAxes, fontsize=fs)
    ax[-1].text( -0.14, 1.05, 'D', transform=ax[-1].transAxes, fontsize=fs)


    fig.savefig(PATH_FIG.joinpath("fig2.pdf"), bbox_inches='tight')



def entropy_effect(base=PATH_BASE, e=8, it=13, dS=10, cut=2, ft=12, fig='', ax='', cax='', km=225):
    # Note that the angles are given in the old format
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]
    ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"{i:03d}", "energy.out")) for i in range(km)])[:,:,0]
    ener = update_energy(ener, ang, e)[:,12:]

    if isinstance(ax, str):
        fig, ax = plt.subplots(1,2, figsize=(10,5))
        fig.subplots_adjust(wspace=0.3)
    ttls = ["Loose Discrimination", "Strict Discrimination"]

    Ks = 10000
    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ent = np.array([-2.0 * np.log(Kmax/Ks)] * ener.shape[1]).T + 1.5 + dS

    cfill = sns.color_palette()[1]
    cdisc = tuple([0.3]*3)
    cnodisc = tuple([0.8]*3)
    al = 0.5

    for i in range(2):
        F = ener + ent
        xlo, xhi = F[:,1:].min(), F[:,1:].max()
        ylo, yhi = F[:,:-1].min(), F[:,:-1].max()

        if i == 0:
            Xlo = np.linspace(xlo, 0, 1000)
            Ylo = Xlo + cut
            Yhi = [yhi] * Xlo.size
            ax[i].fill_between(Xlo, Ylo, Yhi, color=cfill, alpha=al)

            Xhi = np.linspace(xlo, xhi, 1000)
            Ylo = [ylo] * Xhi.size
            Yhi = np.min([(Xhi - cut), np.zeros(Xhi.size)], axis=0)
            ax[i].fill_between(Xhi, Ylo, Yhi, color=cfill, alpha=al)
        else:
            Xlo = np.linspace(xlo, 0, 1000)
            Ylo = np.max([np.zeros(Xlo.size), Xlo + cut], axis=0)
            Yhi = [yhi] * Xlo.size
            ax[i].fill_between(Xlo, Ylo, Yhi, color=cfill, alpha=al)

            Xhi = np.linspace(0, xhi, 1000)
            Ylo = [ylo] * Xhi.size
            Yhi = np.min([(Xhi - cut), np.zeros(Xhi.size)], axis=0)
            ax[i].fill_between(Xhi, Ylo, Yhi, color=cfill, alpha=al)

        X = F[:,1:].ravel()
        Y = F[:,:-1].ravel()

        if i == 0:
            idx = ((Y < 0) & (X - Y >= 2) | (X < 0) & (Y - X >= 2))
        else:
            idx = ((X < 0) & (Y > 0) & (Y - X >= 2) | (X > 0) & (Y < 0) & (X - Y >= 2))
        ax[i].plot(X[idx==False], Y[idx==False], 'o', color=cnodisc, ms=3, label='No Discrimination')
        ax[i].plot(X[idx], Y[idx], 'o', c=cdisc, ms=4, label='Discrimination')

        ax[i].plot([xlo, xhi], [0]*2, '-k')
        ax[i].set_xlim(xlo, xhi)

        ax[i].plot([0]*2, [ylo, yhi], '-k')
        ax[i].set_ylim(ylo, yhi)

        ax[i].set_xlabel(r"$\Delta G_S$ (kT)", fontsize=ft)
        ax[i].set_ylabel(r"$\Delta G_L$ (kT)", fontsize=ft)


    for i in [0,1]:
        ax[i].set_title(ttls[i], fontsize=ft, loc='left')
#       ax[i].annotate(ttls[i], xy=(-12, -1.0+i*2), xycoords='data', xytext=(-10,4.0), horizontalalignment='left',
#                      arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=ft)
#       ax[i].annotate("", xy=(5, -2.0), xycoords='data', xytext=(-3,3.6),
#                      arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1})

    handles = [Line2D([], [], marker='o', color=c, ms=6, linestyle='') for c in [cnodisc, cdisc]]
    lbls = ['No Discrimination', 'Discrimination']
    ax[0].legend(handles, lbls, loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.45, 1.30))

#   cax = inset_axes(ax[0], width="100%", height="10%", loc='upper right',
#                    bbox_to_anchor=(0.0, 1.25, 1.0, 1.0), bbox_transform=ax[0].transAxes, borderpad=0)
#   cbar = fig.colorbar(im, ax=cax, fraction=0.07, pad=0.00, location='top')

#   cax = inset_axes(ax[0], 1, 0.1, loc='upper center')
#   cbar = fig.colorbar(im, ax=cax, fraction=0.046, pad=0.04, location='top')
#   cax.set_axis_off()

    if isinstance(cax, str):
        pass
#       cbar = fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
#       cbar.ax.get_xaxis().labelpad = 12
#       cbar.ax.set_xlabel(r"$\Delta \Delta G$ / kT", rotation=0, fontsize=ft)
    else:
        pass
#       cbar = fig.colorbar(im, cax, orientation='horizontal')
#       cbar.ax.get_xaxis().labelpad = 12
#       cax.set_xlabel(r"$\Delta \Delta G$ / kT", rotation=0, fontsize=ft)
#       cax.xaxis.set_ticks_position('top')
#       cax.xaxis.set_label_position('top')

    fig.savefig(PATH_FIG.joinpath("ent_effect.pdf"), bbox_inches='tight')


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def get_contour(a, xgrid, ygrid):
    b = np.abs(np.diff(a[:,:-1], axis=0)) + np.abs(np.diff(a[:-1,:], axis=1))
    i, j = np.where(b.T)
    X, Y = xgrid[i], ygrid[j]
    idx = np.argsort(X)
    return X[idx], Y[idx]


def plot_rec_space_noseq(dS=20, i0=3, jmax=225, cmap=Oslo_20, lc=Paired_12.hex_colors[10], fig='', ax='', cut=2, ft=12, cax=''):
    if isinstance(ax, str):
        fig = plt.figure(figsize=(12,8))
        gs = GridSpec(3,3, height_ratios=[1,1,1], width_ratios=[2,0.4,1])
        ax = [fig.add_subplot(gs[:,0])] + \
             [fig.add_subplot(gs[i,2]) for i in [0,1,2]]

    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter00", f"Template", "inputs", "ligands.dat"))[:,3][12:]
    Earr = np.arange(0.25, 25.25, 0.25)
    Karr = np.round(10.**np.linspace(1, 5, 300), 2)[:jmax]
    Ks = 10000

    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    ener = np.array([np.load(base.joinpath(f"{i:03d}", "energy.dat.npy")) for i in range(100)])[:,:,:,0]
    ener = np.array([update_energy(ener[i,:jmax,12:], ang, e) for i, e in enumerate(Earr)])
    gap = np.abs(ener[:,:,1:] - ener[:,:,:-1])

    cbar_max = 8

#   return linregress(Karr[gap[:,:,i0].argmax(axis=1)], Earr)[0]

    # Plot ddG as function of K and E
    im = ax[0].contourf(Karr, Earr, gap[:,:,i0], cmap=cmap.mpl_colormap.reversed())

    # Plot optimal line
    xopt, yopt = Karr[gap[:,:,i0].argmax(axis=1)], Earr
    intercept, gradient = linregress(xopt, yopt)[:2]
#   print(gradient)
#   print(linregress(yopt/0.3, gap[:,:,i0].max(axis=1) / (5/180*np.pi)))
    yopt2 = intercept + gradient * xopt
    ax[0].plot(xopt, yopt2, '--', c=sns.color_palette()[1], lw=2.0)



    ax[0].set_xlabel(r"Spring constant, $K$ (kT/nm$^2$)", fontsize=ft)
    ax[0].set_ylabel(r"Energy scale, $\epsilon$ (kT)", fontsize=ft)
    ax[0].set_xlim(10, 10**4)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylim(2, 25)
    if isinstance(cax, str):
        cbar = fig.colorbar(im, ax=ax[0], fraction=0.043, pad=0.04, ticks=np.linspace(0,cbar_max,5))
        cbar.set_label(r"$\Delta \Delta G$ (kT)", rotation=270, fontsize=ft)
        cbar.ax.get_yaxis().labelpad = 20
    else:
#       cbar = fig.colorbar(im, cax, orientation='horizontal')
        cbar = fig.colorbar(im, cax[0], orientation='horizontal', ticks=np.linspace(0,cbar_max,5))
        cbar.ax.get_xaxis().labelpad = 6
        cax[0].set_xlabel(r"$\Delta \Delta G$ (kT)", rotation=0, fontsize=ft)
        cax[0].xaxis.set_ticks_position('top')
        cax[0].xaxis.set_label_position('top')

### We can also link binding entropy to K at this point
    for i, dS in enumerate([-10, 10, 30]):
        ENT = -2.0 * np.log(Karr / Ks) + 1.5 + dS
        ftot = ener[:,:,:] + ENT.reshape(ENT.size, 1)
        rec1 = ((gap[:,:,i0] > cut) & (ftot[:,:,i0+1] < 0))
        rec2 = (gap[:,:,i0] > cut) & (ftot[:,:,i0+1] < 0) & (ftot[:,:,i0] > 0)

        Z = rec1.astype(int) + rec2.astype(int)
        im = ax[i+1].imshow(Z, aspect=Karr.size/Earr.size/1.8, vmin=0, vmax=2, cmap=discrete_cmap(3, cmap.mpl_colormap.reversed()))

        ax[i+1].invert_yaxis()
        ax[i+1].set_xlabel(r"$K$ (kT/nm$^2$)", fontsize=ft)
        ax[i+1].set_ylabel(r"$\epsilon$ (kT)", fontsize=ft)
        ax[i+1].set_xticks(range(Karr.size)[::10])
        ax[i+1].set_xticklabels(np.round(Karr, 0).astype(int)[::10], rotation=90, fontsize=ft)
        ax[i+1].set_yticks(range(Earr.size)[19::20])
        ax[i+1].set_yticklabels(Earr.astype(int)[19::20], fontsize=ft)

        ax[i+1].text(0.04, 0.15, r"$-\Delta S_0={0}$".format(dS), transform=ax[i+1].transAxes, fontsize=ft+1, color='k')

        if i == 1:
            if isinstance(cax, str):
                cax = fig.add_axes([0.9, 0.3, 0.1, 0.5])
                cbar = fig.colorbar(im, ax=cax, fraction=0.043, pad=0.04,
                                    cmap=discrete_cmap(3, cmap.mpl_colormap.reversed()), ticks=np.arange(3)*2/3 + 1/3)
                cbar.ax.set_yticklabels(["No\nDiscrimination", "Loose\nDiscrimination", "Strict\nDiscrimination"], fontsize=ft)
                cbar.ax.get_yaxis().labelpad = 25
            else:
                cbar = fig.colorbar(im, cax[1], cmap=discrete_cmap(3, cmap.mpl_colormap.reversed()), ticks=np.arange(3)*2/3 + 1/3)
                cbar.ax.set_yticklabels(["No\nDisc.", "Loose\nDisc.", "Strict\nDisc."], fontsize=ft)
                cbar.ax.get_yaxis().labelpad = 25
#               cax[1].yaxis.set_ticks_position('right')
#               cax[1].yaxis.set_label_position('right')

        ax[i+1].set_xticks(np.linspace(*ax[i+1].get_xlim(), 4))
        ax[i+1].set_xticklabels([r"$10^{0}$".format(i) for i in range(1, 5)], rotation=0)


def get_alt_energy(cprot, adj, clig, i, k, e, bidx=[1,6,2]):
    return np.array([local_utils.get_energy_simple(cprot[j,i], clig[i].reshape(2,3).T, adj, bidx, k, e) for j in range(len(cprot))])


#ef energy_profile(cprot, karr, clig, c0, adj, fig='', ax=''):
def energy_profile(fig='', ax=''):
    cprot = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", "AllE", "023", f"{i:03d}", "config.out")) for i in range(300)]).reshape(300,35,13,2)
    karr = np.round(10.**np.linspace(1, 5, 300), 2)
    clig = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", "E08", "init_iter12", "000", "inputs", "ligands.dat"))[:,4:]
    c0 = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", "E08", "init_iter12", "000", "inputs", "prot_xy.dat"))
    adj = local_utils.break_links(local_utils.get_adj_from_coord(c0))

    if isinstance(ax, str):
        fig, ax = plt.subplots(3,1)

    i = 15
    e = 10
    col = sns.color_palette()
    pat = '-:'
    for i, i0 in enumerate([15, 17]):
        for j, j0 in enumerate([40, 60, 80, 100]):
            ener = get_alt_energy(cprot, adj, clig, i0, karr[j0], e)
            etot = np.min([ener[:,2], np.ones(300)*(-e)], axis=0)
            edef  = ener[:,0]
            echem = ener[:,1]
            ax[0].plot(karr, edef, pat[i], c=col[j], label=f"lig{i+1}_K={round(karr[j0],0)}")
            ax[1].plot(karr, echem, pat[i], c=col[j])
            ax[2].plot(karr, etot, pat[i], c=col[j])
    for a in ax:
        a.set_xscale('log')
    ax[0].set_ylabel(r"$E_{def}$ given r")
    ax[1].set_ylabel(r"$E_{chem}$ given r")
    ax[2].set_ylabel(r"$E_{tot}$ given r")
    ax[2].set_xlabel(r"value of $K$ used to set $r$")
    ax[0].legend(loc='best', frameon=False)

    return cprot, karr, clig, c0, adj


def get_smooth_envelope(X, Y):
    xgrid = np.linspace(X.min(), X.max(), 100)
    dx = xgrid[1] - xgrid[0]
    ymin, ymax = [], []
    for i, x in enumerate(xgrid):
        y_bin = Y[(X >= x - dx) & (X <= x + dx)]
        ymin.append(y_bin.min())
        ymax.append(y_bin.max())
    return xgrid, np.array(ymin), np.array(ymax)


def specificity_vs_affinity(fig='', ax='', ft=14, i0=3):
    if isinstance(ax, str):
        fig, ax = plt.subplots()

    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    jmax = 225
    Earr = np.arange(0.25, 25.25, 0.25)

    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", "E08", "init_iter01", f"Template", "inputs", "ligands.dat"))[:,3][12:]
    ener = np.array([np.load(base.joinpath(f"{i:03d}", "energy.dat.npy")) for i in range(100)])[:,:jmax,12:,0]
    ener = np.array([update_energy(ener[i], ang, e) for i, e in enumerate(Earr)])
    gap = np.abs(ener[:,:,1:] - ener[:,:,:-1])[:,:,i0]
    F = ener[:,:,i0+1]


### We can also link binding entropy to K at this point
    Ks = 10000
    Karr = np.round(10.**np.linspace(1, 5, 300), 2)[:jmax]
    dS_list = [0, 30]
    for i, dS in enumerate(dS_list):
        ENT = -2.0 * np.log(Karr / Ks) + 1.5 + dS
        ftot = F + ENT.reshape(1, ENT.size)
        ax.scatter(ftot.ravel(), gap.ravel(), alpha=0.03)
#       ax.fill_between(*get_smooth_envelope((ftot).ravel(), gap.ravel()), alpha=0.5, label=r"$-\Delta S = {0}$".format(dS))
        print(dS, "\n", linregress((ftot).ravel(), gap.ravel()))

    xhi, xlo = ax.get_xlim()
    yhi, ylo = ax.get_ylim()

    ax.plot([xlo, xhi], [0,0], '-k')
    ax.plot([0,0], [ylo, yhi], '-k')

    arrow_col = tuple([0.1]*3)
    
    ax.annotate("", xy=(-66, 1.0), xycoords='data', xytext=(-10, 1.0), horizontalalignment='left',
                arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1, "color":arrow_col}, fontsize=ft)
    ax.annotate("", xy=(3,  8.0), xycoords='data', xytext=(3, 0.7), horizontalalignment='right',
                arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1, "color":arrow_col}, fontsize=ft)

    ax.text(-65, 2.0, "High\nAffinity", fontsize=ft)
    ax.text(-28, 2.0, "Low\nAffinity", fontsize=ft)
    ax.text( 7, 1.0, "Low\nSpecificity", fontsize=ft, rotation=0)
    ax.text( 7, 6.0, "High\nSpecificity", fontsize=ft, rotation=0)
    handles = [Line2D([], [], ls='', marker='o', color=sns.color_palette()[i]) for i in range(2)]
    lbls = [r"$-\Delta S_0 = {0}$".format(dS) for dS in dS_list]
    ax.legend(handles, lbls, loc='center right', frameon=False)
    ax.set_xlabel(r"$\Delta G_C$ (kT)", fontsize=ft)
    ax.set_ylabel(r"$\Delta \Delta G$ (kT)", fontsize=ft)
    ax.set_ylim(0, 8.5)


########################################################
### FIG 3
###


def plot_best_df(df, idx, cmap=Vik_20.mpl_colormap, vmin=1, vmax=2, fig='', ax='',
                 bidx=[1,6,2], ec='', cbar=True, side=[0,1], lig=True, cmod=False):

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

        for x, c in zip(xy[bidx], col2):
            c = (c - vmin) / (vmax - vmin)
            if cmod:
                c = c + 0.2 if c < 0.2 else c - 0.2
            ax[i].plot([x[0]], [x[1]], 'o', c=cmap(c), fillstyle='right', ms=13, mec='grey')

        for j, (x, c) in enumerate(zip(xy, col1)):
            c = (c - vmin) / (vmax - vmin)
            if cmod:
                c = c + 0.2 if c < 0.2 else c - 0.2
            ax[i].plot([x[0]], [x[1]], 'o', c=cmap(c), fillstyle='left', ms=13, mec='grey')

            if not isinstance(ec, str):
#               ax[i].quiver([ec0[j,0]], [ec0[j,1]], [x[0] - ec0[j,0]], [x[1] - ec0[j,1]], angles='xy', scale_units='xy', scale=1)
                ax[i].plot([x[0], ec0[j,0]], [x[1], ec0[j,1]], '-k')

        if lig:
            clig = df.loc[idx, f"clig{i+1}"].iloc[0]
            ax[i].plot(*clig.T, 'sk', ms=15, fillstyle='none', mew=2)

    if cbar:
        fig.colorbar(sc, ax=ax[1], fraction=0.046, pad=0.04)


def get_seq_covariance(seq):
#   cov = np.cov(seq.T)
#   np.fill_diagonal(cov, 0)
#   return cov
    N = seq.shape[1]
    cov = np.zeros((N,N), float)
    for i in range(N-1):
        for j in range(i+1, N):
            r = pearsonr(seq[:,i], seq[:,j])[0]
            cov[i,j] = r
            cov[j,i] = r
    return cov


def plot_covariance(df, idx, cmap=Vik_20.mpl_colormap, fig='', ax='', ft=14):
    if isinstance(fig, str):
        fig, ax = plt.subplots()

    cov = get_seq_covariance(np.array(list(df.loc[idx, 'pseq'])))
    cov[np.isnan(cov)] = 0
    vmax = np.max(np.abs(cov))
    vmin = -vmax

    im = ax.imshow(cov, cmap=cmap, vmin=vmin, vmax=vmax)
#   ax.invert_yaxis()
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel("Correlation", rotation=270, fontsize=ft)


### Requires the output of "post-processing.process_test_01",
### using the results of the main FORTRAN program.
### It should be possible to reproduce this within a few hours.
def fig3(df2, ID='w2h1', cmap=Vik_6, ft=10):
    fig = plt.figure(figsize=(8,3))
    fig.subplots_adjust(wspace=2.5)
    gs = GridSpec(1,8)
    ax = [fig.add_subplot(gs[0,i*3+i:(i+1)*3+i]) for i in [0,1]]

    df = df2.loc[(df2.ID==ID)&(df2.ang2==0.785398)]
    idx_list = [df.index, df.loc[(df.gap>3.60)].index]

    for i, idx in enumerate(idx_list):
        plot_covariance(df, idx, fig=fig, ax=ax[i], ft=ft)
        
    xt = np.array([0, 4, 9, 13, 16]) - 0.5
    xl = [f"\n{s}" for s in ["Row 1", "Row 2", "Row 3", "Bind", ""]]

    for a in ax:
        a.set_xticks(xt)
        a.set_yticks(xt)
        lx = a.set_xticklabels(xl, rotation=0, horizontalalignment='left', fontsize=9)
        ly = a.set_yticklabels(xl, rotation=0, verticalalignment='top')

    fig.savefig(PATH_FIG.joinpath(f"fig3.pdf"), bbox_inches='tight', transparent=True)


def best_gap_for_dtheta(i0=17):
    kmax = np.load("/data/Jmcbride/RecProt/Data/NoSeq/K.npy")
    ener = np.load("/data/Jmcbride/RecProt/Data/NoSeq/energy_summary.npy")

    i1 = np.arange(i0-1, max(-1, i0-9), -1)
    gap = np.array([np.abs(ener[:,:,:,i0] - ener[:,:,:,i]).max() for i in i1])
    dtheta = 2.5 * np.arange(gap.size) + 2.5
    return dtheta, gap


def max_gap_task_size(fig='', ax=''):
    if isinstance(fig, str):
        fig, ax = plt.subplots(1,3)
    e_list = [np.load(f"/data/Jmcbride/RecProt/Data/NewCode/free_energy_{i}.npy") for i in range(7)]
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']

    for i, i0 in enumerate([19, 21, 23]):
        i1 = np.arange(i0-1, max(-1, i0-9), -1)
        ang = np.arange(1, i1.size+1) * 2.5
        ax[i].plot(*best_gap_for_dtheta(i0), '-k', label='Ceiling')
        for ID, e in zip(id_list, e_list):
            g = np.array([np.max(np.abs(e[:,i0] - e[:,i])) for i in i1])
            ax[i].plot(ang, g, label=ID)
        ax[i].legend(loc='best', frameon=False)



########################################################
### FIG 4
###


def fig4(out):
    fig, ax = plt.subplots(2,3, figsize=(6,6))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    id_list = ['W1H1', 'W2H1', 'W2H2-', 'W2H2', 'W2H3-', 'W3H1', 'W2H3']
    col = np.array(sns.color_palette())[[0,1,3,2,4]]
    size = np.array([7, 13, 16, 18, 18, 19, 22])
    angles = 2 * (90 - np.arange(30, 95, 5))
    lbls1 = ["Easy Task", "Medium Task", "Hard Task"]
    for i, i0 in enumerate([0, 1, 2]):
        for j, (frac, trav, neigh, nc) in enumerate(out[i0]):
            ms = 8 if j == 3 else 4
            try:
                ax[0,i].plot(angles[1:-1], trav[0][1:-1], '-o', label=id_list[j], c=col[j], ms=ms, mec='k', mew=0.5)
                ax[1,i].plot(angles[1:-1], neigh[0][1:-1], '-o', label=id_list[j], c=col[j], ms=ms, mec='k', mew=0.5)
            except:
                ax[0,i].plot(angles[1:-1], trav[1:-1], '-o', label=id_list[j], c=col[j], ms=ms, mec='k', mew=0.5)
                ax[1,i].plot(angles[1:-1], neigh[1:-1], '-o', label=id_list[j], c=col[j], ms=ms, mec='k', mew=0.5)
        ax[0,i].set_title(lbls1[i])
        ax[1,i].set_title(lbls1[i])


    ymax = [0.95, 0.65]
    for i in [0,1]:
        for j in [0,1,2]:
            ax[i,j].set_xlabel(r"$\theta_C$")
            ax[i,j].set_xticks(np.arange(20, 110, 20))
            ax[i,j].set_xlim(20, 105)
            ax[i,j].set_ylim(-0.03, ymax[i])
    ax[0,0].set_ylabel("Evolvability")
    ax[1,0].set_ylabel("Robustness")


    for a in ax.ravel():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.xaxis.set_minor_locator(MultipleLocator(10))
    
    ax[0,0].legend(loc='upper right', bbox_to_anchor=(4.0, 1.35), frameon=False, ncol=5)

    fs = 16
    for i, b in enumerate('CD'):
        ax[i,0].text( -0.45, 1.05, b, transform=ax[i,0].transAxes, fontsize=fs)

    fig.savefig(PATH_FIG.joinpath(f"fig4.pdf"), bbox_inches='tight')



########################################################
### FIG 5
###


def fig5():
    fig = plt.figure(figsize=(14, 7.0))
    gs = GridSpec(11,3)#2, height_ratios=[.2, 1, 1])
    ax = [fig.add_subplot(gs[1:,0])] + \
         [fig.add_subplot(gs[(i-1)*5+1:i*5+1,1]) for i in [1,2]] + \
         [fig.add_subplot(gs[(i-1)*5+1:i*5+1,2]) for i in [1,2]]
#        [fig.add_subplot(gs[(i-1)*5+i:i*5,2]) for i in [1,2]]
    cax = [fig.add_subplot(gs[0,i]) for i in [0,1]]
    fig.subplots_adjust(hspace=0.05, wspace=0.4)

    spec_aff(fig=fig, ax=ax[:3], cax=cax)
    optimal_mismatch_scaling(fig=fig, ax=ax[3:])

    ax[3].set_xticklabels([''] * len(ax[3].get_xticklabels()))
    for a in ax[3:]:
        a.tick_params(axis='x', direction='in')
    ax[4].invert_yaxis()

#   fs = 16
#   for i, b in zip([0,1,3], 'ABC'):
#       ax[i].text( -0.10, 1.05, b, transform=ax[i].transAxes, fontsize=fs)
    fig.savefig(PATH_FIG.joinpath(f"fig5.pdf"), bbox_inches='tight')


def plot_optimal_line(ax, X, Y, cut, col, offset=0, p='--', d=2):
    idx = range(np.where(Y==cut)[0][0] + 10)
    X, Y = X[idx], Y[idx]
    Yp = np.poly1d(np.polyfit(X, Y, d))(X)
    ax.plot(X, Yp + offset, p, color=col)
    

def spec_aff(base=PATH_BASE, ds=0.0, e=8, it=12, fig='', ax='', cax='', ft=14, km=226, ks=10000):
    if isinstance(ax, str):
        fig = plt.figure(figsize=(7, 4.8))
        gs = GridSpec(3,2, height_ratios=[.2, 1, 1])
        ax = [fig.add_subplot(gs[1:,0])] + \
             [fig.add_subplot(gs[i,1]) for i in [1,2]]
        cax = [fig.add_subplot(gs[0,i]) for i in [0,1]]
        fig.subplots_adjust(hspace=0.05, wspace=0.4)

#   cmap = Oslo_20.mpl_colormap.reversed()
    cmap1 = YlOrBr_9.mpl_colormap
    cmap2 = YlGn_9.mpl_colormap    
    ncont = 16

    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]

    Earr= np.arange(0.25, 25.25, 0.25)
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    i = np.argmin(np.abs(Earr - e))
    ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
    ener = update_energy(ener[:km,:], ang, Earr[i])#/ e

    gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
    gap2 = np.abs(ener[:,14:24] - ener[:,12:22])
    ener = ener[:,12:]
#   ang = ang[12:]

    ent = np.array([1.5 - 2 * np.log(Kmax/ks)] * ener.shape[1]).T + ds
#   ent = np.ones(ent.shape) * -5.5
    fener = np.abs(np.min([ener + ent, np.zeros(ener.shape)], axis=0))
    im = []
    im.append(ax[0].contourf(fener.T[:,::-1], cmap=cmap1, levels=ncont))

    im.append(ax[1].contourf(gap1.T[:,::-1], cmap=cmap2, vmin=0, vmax=e, levels=ncont))
    im.append(ax[2].contourf(gap2.T[:,::-1], cmap=cmap2, vmin=0, vmax=e, levels=ncont))

    # Plot line of optimal mismatch for given flexbility
    ce1 = 'k'
    X = np.arange(km)
    Y = np.argmax(gap1, axis=1)[::-1]
    plot_optimal_line(ax[0], X, Y, 9, ce1, offset=11)
    plot_optimal_line(ax[1], X, Y, 9, ce1)

    Y = np.argmax(gap2, axis=1)[::-1]
    plot_optimal_line(ax[0], X, Y, 0, ce1, offset=2)
    plot_optimal_line(ax[2], X, Y, 0, ce1)

    # Plot line of optimal mismatch for given flexbility for higher energy (e=16)
#   i = 79
#   ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
#   ener = update_energy(ener[:km,:], ang, Earr[i])#/ e
#   gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
#   gap2 = np.abs(ener[:,14:24] - ener[:,12:22])

#   ce2 = 'purple'
#   Y = np.argmax(gap1, axis=1)[::-1]
#   plot_optimal_line(ax[1], X, Y, 9, ce2)

#   Y = np.argmax(gap2, axis=1)[::-1]
#   plot_optimal_line(ax[2], X, Y, 0, ce2)

    # Add colorbar
    clbl = [r"$\Delta G_C$ (kT)", r"$\Delta \Delta G$ (kT)"]
    cticks = [np.arange(0, 30, 8), np.arange(0, 8, 2)]
    for i, j in enumerate([0,2]):
        cbar = fig.colorbar(im[j], cax[i], pad=0.2, orientation='horizontal', ticks=cticks[i])
        cbar.set_label(clbl[i], rotation=0, fontsize=ft-4)
        cbar.ax.get_yaxis().labelpad = 20
        cax[i].xaxis.set_ticks_position('top')
        cax[i].xaxis.set_label_position('top')
    cax[0].set_xticklabels([str(-1 * int(x.get_text())) for x in cax[0].get_xticklabels()])
        
    # Add labels, ticks, etc.
    ylbl = [r"Shape mismatch, $\Delta \theta_{0} ~(^{{\circ}})$",
            r"$\Delta \theta_{0} ~(^{{\circ}})$", r"$\Delta \theta_{0} ~(^{{\circ}})$"]
    xticks = [np.argmin(np.abs(Kmax - 10.**x)) for x in range(1, 5)]
    ylim = [22, 9, 9]
    for i in range(3):
        ax[i].set_xlabel(r"Flexibility, $-\log_{10} K$")
        ax[i].set_ylabel(ylbl[i])
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(range(-4,0,1))
        ax[i].set_ylim(0, ylim[i])
    ax[1].set_xticks([])
    ax[1].set_xlabel("")
    ax[0].set_yticks(np.arange(3, 22, 4))
    ax[0].set_yticklabels(np.arange(-40, 50, 20))
    ax[1].set_yticks(np.arange(0, 12, 4))
    ax[1].set_yticklabels(np.arange(0, 50, 20))
    ax[1].set_yticklabels(["", "20", "40"])
    ax[2].set_yticks(np.arange(1, 12, 4))
    ax[2].set_yticklabels(["-40", "-20", ""])

    fig.savefig(PATH_FIG.joinpath(f"fig6b.pdf"), bbox_inches='tight')


def fit_fn(x, a, b):
    return a + b * x


def optimal_mismatch_scaling(km=226, it=12, sigma=0.3, di=1, fig='', ax=''):
    if isinstance(ax, str):
        fig, ax = plt.subplots(2,1, figsize=(7,5))
        fig.subplots_adjust(hspace=0.3)

    Karr = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    R0 = np.sin(ang-np.pi/3)[12:]
    Earr= np.arange(0.25, 25.25, 0.25)
#   cols = sns.color_palette('husl', n_colors=6)
    cols = sns.color_palette()
    sign = [1, -1]
    e_vals = [2, 20]
    pat = '<>^'
    ms = 8
    al = 0.6

    count = 0
    for i, e in enumerate(Earr):
#       if not ((e == int(e)) and (e/6 == e//6)):
#           continue
        if (int(e) not in e_vals) or (e != int(e)):
            continue
        print(i,e)
        ener0 = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:km,:,0]
        j = np.argmin(np.abs(Earr - e))
        ener = update_energy(ener0[:km,:], ang, Earr[j])[:,12:]

#       gap1 = np.abs(ener[:,11:-di] - ener[:,11+di:])
#       gap2 = np.abs(ener[:,di:12] - ener[:,:12-di])
        gap1 = np.abs(ener[:,11:-di] - ener[:,11+di:])
        gap2 = np.abs(ener[:,di:12] - ener[:,:12-di])

        rmax1 =  R0[11:-di][gap1.argmax(axis=1)]
        rmax2 = -R0[di:12][gap2.argmax(axis=1)]

        for k, (gap, rmax) in enumerate(zip([gap1, gap2], [rmax1, rmax2])):
            imin = np.where(rmax == rmax.max())[0][-1]
            X = e / Karr[imin:] / sigma ** 2
            popt, pcov = curve_fit(fit_fn, X, rmax[imin:] / sigma)
#           print(popt)

            if k == 0:
                lbl = r"$\epsilon={0}$; $\dfrac{{r^*_0}}{{\sigma}} = {1} + {2}\dfrac{{\epsilon}}{{K \sigma^2}}$".format(int(round(e)), round(popt[0],1), round(popt[1],1))
            else:
                lbl = r"$\epsilon={0}$; $\dfrac{{r^*_0}}{{\sigma}} = -{1} - {2}\dfrac{{\epsilon}}{{K \sigma^2}}$".format(int(round(e)), round(popt[0],1), round(popt[1],1))
            ax[k].plot(X, sign[k] * rmax[imin:] / sigma, pat[count], fillstyle='none', ms=ms, alpha=al, c=cols[count], label=lbl)
            ax[k].plot(X, sign[k] * fit_fn(X, *popt), '-', c=cols[count])
        count += 1

    for i, a in enumerate(ax):
        a.set_xlabel(r"$\epsilon / K$")
        a.set_ylabel(r"$\frac{r^*_0}{\sigma}$",rotation=0, fontsize=14, labelpad=15)
        a.set_xlim(0, 0.41)
        a.set_ylim(0, sign[i] * 1.6)
        a.set_yticks(np.arange(0, sign[i] * 2, sign[i] * 0.4))
    fs = 10
    ax[0].legend(loc='lower right',  frameon=False, ncol=1, fontsize=fs,
                 handletextpad=0, columnspacing=1, bbox_to_anchor=(1.02, 0.00))
    ax[1].legend(loc='upper right',  frameon=False, ncol=1, fontsize=fs,
                 handletextpad=0, columnspacing=1, bbox_to_anchor=(1.02, 1.02))
    ax[0].set_yticklabels([''] + list(np.round(np.arange(0.4, 2, 0.4), 1)))
    

#   fig.savefig(PATH_FIG.joinpath(f"fig8.pdf"), bbox_inches='tight')



########################################################
### FIG 6
###


def fig6():
    fig, ax = plt.subplots(1,2, figsize=(8,3.5))
    fig.subplots_adjust(wspace=0.35)
    fs = 12

    X = np.linspace(0, 10, 1000)
    ax[0].plot(X, X, '-k')
    ax[0].fill_between(X, np.zeros(X.size), X, color='red', alpha=0.3)
    ax[0].fill_between(X, X, X+2, color='orange', alpha=0.3)
    ax[0].fill_between(X, X+2, np.ones(X.size)*10, color='green', alpha=0.3)

    [ax[1].plot(X, c/X, '-', c='grey') for c in [2, 7, 16, 28, 45]]

    for a in ax:
        a.set_xlim(0, 10)
        a.set_ylim(0, 10)
        a.set_xticks([])
        a.set_yticks([])
        for direction in ['right', 'top']:
            a.spines[direction].set_visible(False)

    ax[0].set_xlabel("Task Difficulty", fontsize=fs)
    ax[0].set_ylabel("Degrees of Freedom", fontsize=fs)

    ax[1].set_xlabel("Geometric precision", fontsize=fs)
    ax[1].set_ylabel("Dynamic precision", fontsize=fs)

    theta = 46.5
    ax[0].text(0.10, 0.80, "Evolvable\n& Robust", transform=ax[0].transAxes, fontsize=fs)
    ax[0].text(0.15, 0.25, "Molecular Discrimination", transform=ax[0].transAxes, fontsize=fs, rotation=theta)
    ax[0].text(0.40, 0.20, "Underspecified", transform=ax[0].transAxes, fontsize=fs, rotation=theta)

#   ax[1].annotate("More\nDegrees\nof\nFreedom", xy=(0.5,0.5), xycoords='data', xytext=(8.9, 7.0), horizontalalignment='center',
#                  arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=fs)
    ax[1].annotate("More\nDegrees\nof\nFreedom", xy=(0.5,0.5), xycoords='data', xytext=(8.9, 7.0), horizontalalignment='center',
                   arrowprops={'facecolor':'k', "headlength":10, "headwidth":8, "width":1}, fontsize=fs)
#                  arrowprops={'arrowstyle':'->'}, fontsize=fs)

    for i, b in enumerate('AB'):
        ax[i].text( -0.10, 1.05, b, transform=ax[i].transAxes, fontsize=fs+2)

    fig.savefig(PATH_FIG.joinpath(f"fig6.pdf"), bbox_inches='tight')




#######################################################################################
### Presentation figures


def fig1c_alt(base=PATH_BASE, e=8, it=12, ft=14, km=225):
    fig, ax = plt.subplots(figsize=(8,4))

    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]
    ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E{e:02d}", f"init_iter{it:02d}", f"{i:03d}", "energy.out")) for i in range(km)])[:,:,0]
#   ener = np.array([np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"AllE", f"023", f"{i:03d}", "energy.out")) for i in range(300)])[:,:,0]
    ener = update_energy(ener, ang, e) / e
    gap = ener[:,:-2] - ener[:,2:]

    X = 2 * (90 - ang[:-2]*180/np.pi)
    theta = np.pi/3
    j0 = np.argmin(np.abs(ang - theta))
    ax.plot(Kmax, gap[:,14], '-k', lw=2)


    ax.set_xlim(7, 50000)
    ax.set_xscale('log')
    ax.set_xlabel(r"Spring Constant, $K$ (kT/nm$^{2}$)", fontsize=ft)
#   ax.set_ylabel("Binding Energy gap\n $\\Delta G_B - \\Delta G_A$", fontsize=ft)
    ax.set_ylabel(r"$\Delta \Delta G_{LS}$ ($\epsilon$ kT)", fontsize=ft)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(ft)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ft)

    fig.savefig(PATH_FIG.joinpath("PRES_fig1c.pdf"), bbox_inches='tight')







