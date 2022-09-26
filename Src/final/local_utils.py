import os
from pathlib import Path

import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
from palettable.scientific.diverging import Berlin_20, Roma_20, Berlin_10, Vik_20, Vik_6
from palettable.scientific.sequential import Bilbao_20
import pandas as pd
from scipy.optimize import curve_fit
from scipy.spatial import distance_matrix
import seaborn as sns
from sklearn.cluster import DBSCAN

import new_figs

PATH_BASE = Path("/molahome/jmcbride/RecProt")


def plot_heatmap(Zi, xbins, ybins, xlbl, ylbl, fig='', ax=''):
    if isinstance(ax, str):
        fig, ax = plt.subplots()
    im = ax.imshow(Zi.T)
    ax.invert_yaxis()
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_xticks(range(xbins.size)[::2])
    ax.set_yticks(range(ybins.size)[::2])
    ax.set_xticklabels(np.round(xbins, 2)[::2], rotation=90)
#   ax.set_xticklabels(np.round(np.log10(xbins), 1))
    ax.set_yticklabels(np.round(ybins, 1)[::2])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_combo(bins, combo, totals, recog, typ='dens', istart=0, istep=1, save=False, itst=0):
    fig, ax = plt.subplots()
    xlbl, ylbl = combo
    if typ in ['prob', 'tot']:
        cs = f"{xlbl}_{ylbl}_dens"
    else:
        cs = f"{xlbl}_{ylbl}_{typ}"
    cs0 = f"{xlbl}_{ylbl}_{typ}"
    try:
        if typ == "dens":
            R = np.nansum(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
            plot_heatmap(np.log10(R), bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax)
        elif typ == 'tot':
            T = np.nansum(np.array([x for x in totals[cs][istart::istep] if len(x)]), axis=0)
            plot_heatmap(np.log10(T), bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax) 
        elif typ == 'prob':
            R = np.nansum(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
            T = np.nansum(np.array([x for x in totals[cs][istart::istep] if len(x)]), axis=0)
            plot_heatmap(R / T, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax) 
        elif typ == 'Z':
            R = np.nanmean(np.array([x for x in recog[cs][istart::istep] if len(x)]), axis=0)
            plot_heatmap(R, bins[xlbl], bins[ylbl], xlbl, ylbl, fig, ax) 

        if save:
            for d in ['All', f'Test{itst:02d}', f"{cs}"]:
                p = Path(f"/data/Jmcbride/RecProt/N13/Figures/{d}/")
                p.mkdir(parents=True, exist_ok=True)
                fig.savefig(f"/data/Jmcbride/RecProt/N13/Figures/{d}/e{istart}_{itst:02d}_{cs0}.pdf", bbox_inches='tight')
                fig.savefig(f"/data/Jmcbride/RecProt/N13/Figures/{d}/e{istart}_{itst:02d}_{cs0}.png", bbox_inches='tight')
    except Exception as e:
        print(f"Error at {itst}, {cs}:\n{e}")

    plt.close()


def make_all_plots():
    for i in range(14):
        ts = time.time()
        print(i)
        bins, combos, totals, recog = pickle.load(open(f"../N13/test_{i:02d}_2dplots.pickle", 'rb'))
        for j, c in enumerate(combos):
            for k in range(2):
                plot_combo(bins, c, totals, recog, typ='tot', save=True, itst=i, istart=k, istep=2)
                plot_combo(bins, c, totals, recog, typ='dens', save=True, itst=i, istart=k, istep=2)
                plot_combo(bins, c, totals, recog, typ='prob', save=True, itst=i, istart=k, istep=2)
                if c[1] != 'gap':
                    plot_combo(bins, c, totals, recog, typ='Z', save=True, itst=i, istart=k, istep=2)

        print(f"Time taken: {(time.time()-ts)/60} minutes")


def load_pickle(path):
    return pickle.load(open(path, 'rb'))


def load_binding_energy():
#   base = Path("/molahome/jmcbride/RecProt/Results/BigSmall/Iter01")
#   base = Path("/molahome/jmcbride/RecProt/Results/LocalStruct/Rules11")
    base = Path("/molahome/jmcbride/RecProt/Results/LocalStruct/Rules12")
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
    size = [7, 13, 16, 18, 18, 19, 22]
    e_list = []
    for ID, s in zip(id_list, size):
        if ID in ['w1h1', 'w2h1']:
            e_list.append(np.load(base.joinpath(f"{ID}/293/free_energy.npy")))
        else:
            with Pool(60) as pool:
                imax = int(np.ceil(2.**(s+3)/100000))
                print(ID, s, imax)
                e_list.append(np.concatenate(list(pool.map(np.load, [base.joinpath(f'{ID}/S{i:03d}/free_energy.npy') for i in range(imax)], 1))))
    return e_list


def load_equil_coord():
    base = Path("/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/w2h1/Template")
    return np.loadtxt(base.joinpath("inputs", "prot_xy.dat"))


def load_solutions(base=Path("/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/")):
#   base=Path("../../Results/OpenOpenSize/Iter02/")
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
    df1 = []
    df2 = []
    for i, ID in enumerate(id_list):
        print(i, ID)
        try:
            df1.append(pd.read_pickle(base.joinpath(f"{ID}", "rec1.pickle")))
        except:
            pass
        try:
            df2.append(pd.read_pickle(base.joinpath(f"{ID}", "rec2.pickle")))
        except:
            pass

    df1 = pd.concat(df1, ignore_index=True)
    df2 = pd.concat(df2, ignore_index=True)
    return df1, df2


def get_energy_simple(cprot, clig, adj, bidx, k, e, l0=1, sigma=0.3):
    Edef = 0.0
    for i, j in zip(*np.where(adj)):
        if i < j:
            dr = cprot[i] - cprot[j]
            Edef += k * (dr.dot(dr)**0.5 - l0)**2
#           print(i, j, dr.dot(dr)**0.5)
#           print(k * (dr.dot(dr)**0.5 - l0)**2)

    Echem = 0.0
    for i in range(3):
        dr = cprot[bidx[i]] - clig[i]
        Echem -= e * np.exp(-(dr.dot(dr))/sigma**2)
#       print(i, dr.dot(dr)**0.5)
#       print(e * np.exp(-(dr.dot(dr))**0.5/sigma))
        
#   print(Edef, Echem)
    return [Edef, Echem, Edef + Echem]


def get_rigid_binding_energy(angle, e, a0=np.pi/3, sigma=0.3):
    chord = 2 * np.sin((a0 - angle)/2)
    if angle < a0:
        return -e
    elif abs(angle - a0) < 0.001:
        return -e * 3
    elif angle > a0:
        return np.sum([- e * np.exp(-(r * chord)**2/sigma**2) for r in [0, 0, 1]])


def plot_covariance(seq, cmap=Vik_20.mpl_colormap, fig='', ax=''):
    if isinstance(fig, str):
        fig, ax = plt.subplots()

    cov = np.cov(seq.T)
    np.fill_diagonal(cov, 0)
    vmax = np.max(np.abs(cov))
    vmin = -vmax

    im = ax.imshow(cov, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.invert_yaxis()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def get_cov_df(df, idx):
    seq = np.array(list(df.loc[idx, 'pseq']))
    return get_cov_seq(seq)


def get_cov_seq(seq):
    cov = np.cov(seq.T)
    np.fill_diagonal(cov, 0)
    return cov


def plot_cov_vs_dist(d, c, alg='mean'):
    D = np.round(d, 2)
    C = np.abs(c)
    duniq = np.array(sorted(np.unique(D)))
    if alg == 'mean':
        out = np.array([np.mean(C[D==x]) for x in duniq])
    elif alg == 'max':
        out = np.array([np.max(C[D==x]) for x in duniq])
    return duniq, out


def get_cov_vs_d(df, idx, coord, alg='max', bidx=[1,2,6]):
    D = distance_matrix(coord, coord)
    cov = get_cov_df(df, idx)
    X, Y = [], []
    for i, j in zip(*np.triu_indices(len(coord), 1)):
        if i in bidx and j in bidx:
            continue
        if not i in bidx and not j in bidx:
            continue
        X.append(D[i,j])
        Y.append(cov[i,j])
    return plot_cov_vs_dist(X, Y, alg)


def plot_all_cov_vs_d(df2):
    c_list = [np.loadtxt(f'../../Results/BigSmall/Iter01/{ID}/Template/inputs/prot_xy.dat') for ID in id_list]
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-']

    
    for ID, c in zip(id_list[:5], c_list):
        if ID == 'w1h1':
            bidx = [0,1,3]
        else:
            bidx = [1,2,6]
        plt.plot(*get_cov_vs_d(df2, (df2.ID==ID)&(df2.ang2==0.785398)&(np.abs(df2.gap)>2), c, 'mean', bidx), '-o', label=ID)


def update_energy(ener, ang, seq, emat):
    ener = ener.copy()
    theta = np.pi/3
    seq = seq.copy() - 1

    idx0 = np.zeros(ang.size, bool)

    idx = ang > theta
    idx0 = idx0 + idx
    count = np.sum(idx)
    for i in np.where(idx)[0]:
        chord = 2 * np.sin(ang[i] - theta)
        edx = np.exp(-chord**2/0.3**2)
        emed = emat[seq[1]] + max(emat[seq[0]], emat[seq[1]]) + min(emat[seq[0]], emat[seq[1]]) * edx 
        ener[:,i] = np.min([ener[:,i], np.zeros(ener.shape[0]) - emed], axis=0)

    emin = max(emat[seq[0]], emat[seq[1]])
    idx = ang < theta
    idx0 = idx0 + idx
    count += np.sum(idx)
    ener[:,idx] = np.min([ener[:,idx], np.zeros((ener.shape[0], idx.sum())) - emin], axis=0)

    emax = np.sum(emat[seq])
    idx = np.abs(ang - theta) < 0.001
    idx0 = idx0 + idx
    count += np.sum(idx)
    ener[:,idx] = np.min([ener[:,idx], np.zeros((ener.shape[0], idx.sum())) - emax], axis=0)

    return ener


def load_noseq_ener():
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
    angles = np.arange(2.5, 90, 2.5) * np.pi / 180
    emat = np.array([1.5, 6])
    base = Path("/molahome/jmcbride/RecProt/Results")
    seq = np.loadtxt(base.joinpath('NoSeq/Size/w1h1/Template/inputs/prot_seq.dat'))[:,-3:].astype(int)
    ener = np.array([np.load(base.joinpath(f"NoSeq/Size/{ID}/energy.dat.npy")).reshape(300,8,35,4)[:225,:,:,0] for ID in id_list])

    for i, ID in enumerate(id_list):
        for j, s in enumerate(seq):
            ener[i,:,j,:] = update_energy(ener[i,:,j,:], angles, s, emat)

    return ener


def best_gap_for_dtheta(i0=17):
    kmax = np.load("/data/Jmcbride/RecProt/Data/NoSeq/K.npy")
    ener = np.load("/data/Jmcbride/RecProt/Data/NoSeq/energy_summary.npy")

    i1 = np.arange(i0-1, max(-1, i0-9), -1)
    gap = np.array([np.abs(ener[:,:,:,i0] - ener[:,:,:,i]).max() for i in i1])
    dtheta = 2.5 * np.arange(gap.size) + 2.5
    return dtheta, gap


def plot_gap_vs_dtheta_by_shape(e_list, i0=17):
    fig, ax = plt.subplots()
    i1 = np.arange(i0-1, max(-1, i0-9), -1)
    ang = np.arange(1, i1.size+1) * 2.5
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']

    ax.plot(*best_gap_for_dtheta(i0), '-k', label='Ceiling')
    for ID, e in zip(id_list, e_list):
        g = np.array([np.max(np.abs(e[:,i0] - e[:,i])) for i in i1])
        ax.plot(ang, g, label=ID)
    ax.legend(loc='best', frameon=False)


def plot_d_vs_c_task_shape(s_list, e_list, d_list, i0=19, di=2, cut=2, alg='max', fig='', ax=''):
    if isinstance(fig, str):
        fig, ax = plt.subplots()
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
    col = sns.color_palette()
    for i, ID in enumerate(id_list[:-2]):
        G = e_list[i][:,i0] - e_list[i][:,i0-di]
        idx = (G < -cut) & (e_list[i][:,i0] < 0)
        print(i, ID, np.sum(idx), np.sum(idx) / len(idx))

#       idx = np.array([i for i in np.argsort(G) if idx[i]])[:500]
#       print(G[idx].min(), G[idx].max(), G.min(), G.max())
#       print(len(idx), '\n')

        if np.sum(idx) == 0:
            print(f'No solutions for {ID}')
            continue

        cov = get_cov_seq(s_list[i][idx,:-3])

        idx_ran = np.random.choice(range(len(s_list[i])), size=len(idx))
        cov_ran = get_cov_seq(s_list[i][idx_ran,:-3])

        if i == 0:
            bidx = [0,1,3]
        elif i == 5:
            bidx = [2,3,9]
        else:
            bidx = [1,2,6]

        X, Y, Y2 = [], [], []
        for j, k in zip(*np.triu_indices(s_list[i].shape[1]-3, 1)):
            if j in bidx and k in bidx:
                continue
            if not j in bidx and not k in bidx:
                continue
            X.append(d_list[i][j,k])
            Y.append(cov[j,k])
            Y2.append(cov_ran[j,k])
        D, Cm = plot_cov_vs_dist(X, Y, alg)
        D, Cmr = plot_cov_vs_dist(X, Y2, alg)
        ax.plot(D, Cm - Cmr, label=ID, c=col[i])
#       ax.plot(D, Cm, label=ID, c=col[i])

    ax.plot([1,5], [0,0], '-k')

    ax.legend(loc='best', frameon=False)


def plot_all_above(s_list, e_list, d_list, i0=19, cut=2, alg='max'):
    fig, ax = plt.subplots(1, 3, sharex=True, sharey=True)
    for i, di in enumerate([4,2,1]):
        plot_d_vs_c_task_shape(s_list, e_list, d_list, i0=i0, di=di, cut=cut, alg=alg, fig=fig, ax=ax[i])
#       plot_d_vs_c_task_shape(s_list, e_list, d_list, i0=i0, di=di, cut=di, alg=alg, fig=fig, ax=ax[i])
        

def plot_solution_gap_cum_dist(ener, i=17, j=15, cut=2, dx=0.01):
    gap = ener[:,i] - ener[:,j]
    dG = ener[:,i]
    gap_sol = gap[(gap <= -cut) & (dG<0)]
    if not len(gap_sol):
        return
    edges = np.arange(gap_sol.min(), gap_sol.max()+dx, dx)
    hist, edges = np.histogram(gap_sol, bins=edges)
    hist = hist / hist.sum()
    X = edges[:-1] + 0.5 * np.diff(edges[:2])
    plt.plot(X, np.cumsum(hist))


def update_energy_def(ener, ang, e):
    new_ener = ener.copy()
    theta = np.pi/3

    idx = ang < theta
    new_ener[:,idx,0] = np.min([ener[:,idx,0], np.zeros((ener.shape[0], idx.sum())) - e], axis=0)
    idx2 = ener[:,idx,0] > new_ener[:,idx,0]
    new_ener[:,idx,1][idx2] = 0 

    idx = ang == theta
    new_ener[:,idx,0] = np.min([ener[:,idx,0], np.zeros((ener.shape[0], idx.sum())) - e*3], axis=0)
    idx2 = ener[:,idx,0] > new_ener[:,idx,0]
    new_ener[:,idx,1][idx2] = 0 

    idx = ang > theta
    new_ener[:,idx,0] = np.min([ener[:,idx,0], np.zeros((ener.shape[0], idx.sum())) - e*2], axis=0)
    idx2 = ener[:,idx,0] > new_ener[:,idx,0]
    new_ener[:,idx,1][idx2] = 0 

    return new_ener


def get_adj_from_coord(c, d0=1.01):
    dist = distance_matrix(c, c)
    np.fill_diagonal(dist, dist.max())
    return (dist<=d0).astype(int)


def break_links(adj, idx=[[1,2]]):
    for i, j in idx:
        adj[i,j] = 0 
        adj[j,i] = 0 
    return  adj 


def fit_fn(x, a, b):
    return a + b * x


def testing_tsvi(ax, dS=20, jmax=225, km=226, e=8, it=12, ks=10000, ds=0):
    Karr = np.round(10.**np.linspace(1, 5, 300), 2)[:km]
    ang = np.loadtxt(PATH_BASE.joinpath("Results", "NoSeq", f"E04", f"init_iter{it:02d}", f"Template", "inputs", "ligands.dat"))[:,3]

    Earr= np.arange(0.25, 25.25, 0.25)
    base = PATH_BASE.joinpath("Results", "NoSeq", "AllE")
    i = np.argmin(np.abs(Earr - e))
    ener = np.load(base.joinpath(f"{i:03d}", "energy.dat.npy"))[:,:,0]
    ener = new_figs.update_energy(ener[:km,:], ang, Earr[i])#/ e

#   gap1 = np.abs(ener[:,23:-2] - ener[:,25:])
    gap1 = np.abs(ener[:,23:-1] - ener[:,24:])
#   gap2 = np.abs(ener[:,14:24] - ener[:,12:22])
    gap2 = np.abs(ener[:,13:24] - ener[:,12:23])
    ener = ener[:,12:]
#   ang = ang[12:]

    ent = np.array([1.5 - 2 * np.log(Karr/ks)] * ener.shape[1]).T + ds
#   ent = np.ones(ent.shape) * -5.5
    fener = np.abs(np.min([ener + ent, np.zeros(ener.shape)], axis=0))

    R0 = np.sin(np.pi/3-ang)

#   fig, ax = plt.subplots()
#   rmax = -R0[23:-1][gap1.argmax(axis=1)]
    rmax = R0[13:24][gap2.argmax(axis=1)]

    imin = np.where(rmax == rmax.max())[0][-1]
#   X = -np.log10(Karr[imin:])
#   ax.plot(X, rmax[imin:] / 0.3)
#   ax.set_xlabel(r"$-\log_{10} K$")
    ax.set_xlabel(r"$\epsilon / K$")
    ax.set_ylabel(r"$R_0 / \sigma$")

    X = e / Karr[imin:]
    X2 = -np.log10(Karr[imin:])
    popt, pcov = curve_fit(fit_fn, X, rmax[imin:] / 0.3)
    print(popt)
    plt.plot(X, rmax[imin:] / 0.3)
    ax.plot(X, fit_fn(X, *popt), ':',)

#   popt, pcov = curve_fit(fit_fn, X, rmax[imin:] / 0.3)
#   print(popt)
#   ax.plot(X, fit_fn(X, *popt), ':',)

#   plt.plot(X, np.exp(X))
#   plt.plot(1 + e / Karr- rmax / 0.3)
#   plt.plot((rmax / 0.3 - 1) * Karr / e)
#   plt.plot(np.log(Karr), np.exp(-np.log(Karr))

#   for i0 in range(gap1.shape[1]):
#       print(ang[i0] * 180 / np.pi, R0[i0])
#       print(gap1[:,i0].argmax())
#       print(gap1[:,i0].max())
#       kmax = Karr[gap1[:,i0].argmax()]
#       ax.plot(Karr, gap1[:,i0])

#       ax.plot(1 + epsilon / kmax)
#       popt, pcov = curve_fit(fit_fn, epsilon / Karr, -R0[i0] / 0.3)
#       print(popt[0], '\n')




