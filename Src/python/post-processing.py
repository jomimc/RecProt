"""
Tools for processing the output filess of the main FORTRAN program
"""
import argparse
import glob
from pathlib import Path
import pickle

import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from multiprocessing import Pool
import numpy as np
from palettable.scientific.diverging import Berlin_20, Roma_20, Berlin_10, Vik_20
import pandas as pd
from scipy.signal import argrelmin
import seaborn as sns

import create_inputs
import get_entropy as GE
import plots
import rec_io
from rec_io import PATH_BASE
import rec_utils as utils


N_PROC = 40


#########################################################
### Load initial results

### These came from running 180 ligands that only differ
### by angle. Only one protein spring, and only one
### chemical bond strength.


def load_results(path):
    data = rec_io.load_output(path)
    etot, edef, echem, einit = data.T
    clig = data[:,2:8]
    cprot = data[:,13:]
    seq = data[:,8]
    return clig, cprot, etot, edef, echem, einit, seq


def load_results_new(base):
    etot, edef, echem, einit = np.loadtxt(base.joinpath("energy.out")).T
    clig = rec_io.parse_ligands(base.joinpath("inputs", "ligands.dat"))[2]
    cprot = np.loadtxt(base.joinpath("config.out"))
    seq = np.loadtxt(base.joinpath("inputs", "prot_seq.dat"))
    return clig, cprot, etot, edef, echem, einit, seq


def load_equil_coord(base):
#   return np.loadtxt(PATH_BASE.joinpath("Src", "test_angle", "inputs", "prot_xy.dat"))
    return np.loadtxt(base.joinpath("inputs", "prot_xy.dat"))


### Print out some candidates for further testing

def angle_pairs_equal_deformation():
    clig, cprot, etot, edef, echem, einit, seq = load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    deform = np.arange(0.1, 0.6, 0.1)
    total = [-4.4, -4.2, -4.0, -3.8]
    angles = np.arange(180) * np.pi / 2 / 180
    return etot, edef, echem, angles

    for d in deform:
        idx = argrelmin(np.abs(edef-d))
        print(f"For e_def = {d}, closest match is found for...")
#   for d in total:
#       idx = argrelmin(np.abs(etot-d))
#       print(f"For e_tot = {d}, closest match is found for...")
        print("Deformation energy: ", edef[idx])
        print("Total energy: ", etot[idx])
        print("Angles: ", (angles*180/np.pi)[idx])
        print("")



#########################################################
### Create inputs for tests


def clamp_open_test():
    clig, cprot, etot, edef, echem, einit, seq = load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    deform = np.arange(0.1, 0.6, 0.1)
    total = [-4.4, -4.2, -4.0, -3.8]
    angles = np.arange(180) * np.pi / 2 / 180

    ligands = []
    init_config = []
    ligand_colors = [2, 2, 2]

    for d in deform:
        idx = argrelmin(np.abs(edef-d))
        for i in idx[0]:
            ligands.append(ligand_colors + [angles[i]])
            c = cprot[i]
            c[1::2] = c[1::2] + 0.001
            init_config.append(c)

    for d in total:
        idx = argrelmin(np.abs(etot-d))
        for i in idx[0]:
            ligands.append(ligand_colors + [angles[i]])
            c = cprot[i]
            c[1::2] = c[1::2] + 0.001
            init_config.append(c)

    return ligands, np.array(init_config)


def write_test_ligands_01():
    ligands, cprot = clamp_open_test()

    path = PATH_BASE.joinpath("Results", "ClampOpen", "Template", "inputs", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)

    path = PATH_BASE.joinpath("Results", "ClampOpen", "Template", "inputs",  "prot_xy_init.dat")
    np.savetxt(path, cprot)


def open_open_test():
    clig, cprot, etot, edef, echem, einit, seq = load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    angles_wanted = np.arange(0, 35, 5) * np.pi / 180
    angles = np.arange(180) * np.pi / 2 / 180

    ligands = []
    init_config = []
    ligand_colors = [2, 2, 2]

    for a in angles_wanted:
        i = np.argmin(np.abs(angles-a))
        ligands.append(ligand_colors + [a])
        c = cprot[i]
        c[1::2] = c[1::2] + 0.001
        init_config.append(c)

    return ligands, np.array(init_config)


def write_test_ligands_02():
    ligands, cprot = open_open_test()

    path = PATH_BASE.joinpath("Results", "OpenOpen", "Template", "inputs", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)

    path = PATH_BASE.joinpath("Results", "OpenOpen", "Template", "inputs",  "prot_xy_init.dat")
    np.savetxt(path, cprot)


def sym_asym_test():
    clig, cprot, etot, edef, echem, einit, seq = load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    angles = np.arange(180) * np.pi / 2 / 180
    new_lig = [[1, 1, 2], [1, 2, 1]]
    new_angles = np.arange(0,180,30) * np.pi / 2 / 180
    ligands = []
    init_config = []
    for a in new_angles:
        i = np.argmin(np.abs(angles-a))
        for l in new_lig:
            ligands.append(l + [a])
            c = cprot[i]
            c[1::2] = c[1::2] + 0.001
            init_config.append(c)

    return ligands, np.array(init_config)


def write_test_ligands_03():
    ligands, cprot = sym_asym_test()

    path = PATH_BASE.joinpath("Results", "SymAsym", "Template", "inputs", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)

    path = PATH_BASE.joinpath("Results", "SymAsym", "Template", "inputs",  "prot_xy_init.dat")
    np.savetxt(path, cprot)


def clamp_open_sym_asym_test(switch=0):
    clig, cprot, etot, edef, echem, einit, seq = load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    angles = np.arange(180) * np.pi / 2 / 180
    deform = np.arange(0.1, 0.6, 0.1)
    if switch:
        new_lig = [[1, 1, 2], [1, 2, 1]]
    else:
        new_lig = [[1, 2, 1], [1, 1, 2]]

    ligands = []
    init_config = []
    for d in deform:
        idx = argrelmin(np.abs(edef-d))[0]
        for i, l in zip(idx, new_lig):
            ligands.append(l + [angles[i]])
            c = cprot[i]
            c[1::2] = c[1::2] + 0.001
            init_config.append(c)

    return ligands, np.array(init_config)


def write_test_ligands_04():
    ligands, cprot = clamp_open_sym_asym_test(1)
    path = PATH_BASE.joinpath("Results", "ClampOpenSymAsym", "Template", "inputs", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)

    ligands, cprot = clamp_open_sym_asym_test(0)
    path = PATH_BASE.joinpath("Results", "ClampOpenAsymSym", "Template", "inputs", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)

    path = PATH_BASE.joinpath("Results", "ClampOpenSymAsym", "Template", "inputs",  "prot_xy_init.dat")
    np.savetxt(path, cprot)

    path = PATH_BASE.joinpath("Results", "ClampOpenAsymSym", "Template", "inputs",  "prot_xy_init.dat")
    np.savetxt(path, cprot)




#########################################################
### Analyse test results


def load_test_results(base):
    # Load raw results
    clig, cprot, etot, edef, echem, einit, seq = load_results_new(base)
    params = rec_io.parse_params(base.joinpath("inputs", "parameters.txt"))
    kmat = rec_io.get_spring_const(params)
    emat = rec_io.get_chem_const(params)

    nprot = params['nprot']
    nlig = params['nligand']
    Naa = params['namino']


    # Reshape arrays
    clig = clig[:nlig].reshape(nlig, 3, 2)
    cprot = cprot.reshape(nprot, nlig, Naa, 2)
    etot = etot.reshape(nprot,nlig)
    edef = edef.reshape(nprot,nlig)
    echem = echem.reshape(nprot,nlig)
    einit = einit.reshape(nprot,nlig)

    lseq, angles = rec_io.parse_ligands(base.joinpath("inputs", "ligands.dat"))[:2]

    # Checking that no protein bonds are crossing
    ec = load_equil_coord(base)
    # Don't run "break_links" yet, as this allows to check that the
    # ligand binding site does not have any edges crossed
    adj = utils.get_adj_from_coord(ec)
    is_valid = check_triangle_orientation(ec, cprot, adj)

    # Checking that the solution is better than one of the trivial ones
    enodef = np.array([[no_def_energy(s[Naa:], ls, a, emat) for s in seq] for ls, a in zip(lseq, angles)]).T
   
    # Calculating entropy
    kstiff = np.max(kmat) * 10
    idx_break = [[params["ibind"][0], params["ibind"][2]]]
    # Break the bond at the binding site
    adj = utils.break_links(adj, idx_break)
    ent = np.array([GE.entropy_diff(ec, adj, kmat, s, idx0=params["ibind"], kstiff=kstiff) for s in seq])

    return clig, cprot, etot, edef, echem, einit, lseq, angles, kmat, emat, is_valid, enodef, ent, seq


def plot_gap_01(etot, angles, enodef):
    fig, ax = plt.subplots(2,1)
    for i in range(int(etot.shape[1]/2)):
        lbl = '_'.join([str(round(x,1)) for x in angles[i*2:(i+1)*2]])
        idx = (etot[:,i*2] < enodef[:,i*2]) & (etot[:,i*2+1] < enodef[:,i*2+1])
        gap = etot[idx,i*2] - etot[idx,i*2+1]
        j = i//5
        sns.distplot(gap, label=lbl, ax=ax[j])

    for a in ax:
        a.legend()
        a.set_xlabel("Recognition gap")


def plot_best_01(fig, ax, clig, cprot, seq, idx, ec=''):
    ax.plot(*ec.T, 'sk', fillstyle='none', ms=10)
    if isinstance(ec, str) or 1:
        ec = load_equil_coord(PATH_BASE.joinpath("Src", "test_angle"))
#       ec[:,1] = ec[:,1] - 0.5 * 3**0.5
        ax.plot(*ec.T, 'ok', fillstyle='none', ms=10)
    ax.plot(*clig.reshape(3,2).T, '*', fillstyle='none', ms=20, c='green')

    sc = ax.imshow(np.array([[0.1,0.2],[0.1, 0.2]]), cmap=Berlin_20.mpl_colormap, vmin=1, vmax=2)
    ax.invert_yaxis()
    sc.set_visible(False)

    xy = cprot[idx].mean(axis=0).reshape(13,2)
    col = seq[idx].mean(axis=0)[:13]
    for x, c in zip(xy, col):
        ax.plot([x[0]], [x[1]], 'o', c=Berlin_20.mpl_colormap(c-1), fillstyle='left', ms=10)
#   sc = ax.scatter(*xy.T, c=col, vmin=1, vmax=2, cmap=Berlin_20.mpl_colormap, marker=MarkerStyle('o', fillstyle='left'))

    col = seq[idx].mean(axis=0)[13:]
    for x, c in zip(xy[[1,6,2]], col):
        ax.plot([x[0]], [x[1]], 'o', c=Berlin_20.mpl_colormap(c-1), fillstyle='right', ms=10)
#   sc = ax.scatter(*xy[[0,1,6]].T, c=col, vmin=1, vmax=2, cmap=Berlin_20.mpl_colormap, marker=MarkerStyle('o', fillstyle='right'))
#   fig.colorbar(Berlin_20.mpl_colormap, ax=ax, vmin=1, vmax=2)
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)


def plot_all_best_01(clig, cprot, seq, etot, switch=0):
    ec = np.loadtxt(PATH_BASE.joinpath("Results", "ClampOpen", "Test02", "inputs", "prot_xy_init.dat")).reshape(18,13,2)
    if switch == 0:
        istart, itot = 0, 5
        lo = np.array([-0.015, 0, -0.01, -0.03 , -0.03])
        hi = np.array([0, 0.04, 0.05, 0.06, 0.1])
    else:
        istart, itot = 10, 4
#       lo = np.array([-0.015, -0.02, -0.05, -0.1])
#       hi = np.array([0, 0.04, 0.05, 0.08])
#       lo = np.array([-0.12, -0.19, -0.93, -1.57])
#       hi = np.array([0.35, 0.82, 3.95, 7.85])
        lo = np.array([-0.21, -0.45, -1.20, -2.00])
        hi = np.array([0.18, 0.90, 1.15, 6.05])
#   lo, hi = lo * 10, hi * 10

    fig, ax = plt.subplots(2,itot)
    for i in range(itot):
        j = istart + i*2
        k = istart + i*2 + 1
        plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, (etot[:,j] - etot[:,k]) < lo[i], ec=ec[j])
        print(j, np.sum((etot[:,j] - etot[:,k]) < lo[i]))
        plot_best_01(fig, ax[1,i], clig[k], cprot[:,k], seq, (etot[:,j] - etot[:,k]) > hi[i], ec=ec[k])
        print(k, np.sum((etot[:,j] - etot[:,k]) > hi[i]))


def plot_gap_02(etot, angles, enodef, is_valid):
    fig, ax = plt.subplots()
    for i in range(etot.shape[1]-1):
        lbl = '_'.join([str(round(x,2)) for x in angles[i:i+2]])
        idx = (etot[:,i] < enodef[:,i]) & (etot[:,i+1] < enodef[:,i+1]) & is_valid[:,i] & is_valid[:,i+1]
        gap = etot[idx,i] - etot[idx,i+1]
        sns.distplot(gap, label=lbl)
    ax.legend()
    ax.set_xlabel("Recognition gap")


def plot_all_best_02(clig, cprot, seq, etot, enodef):
    ec = np.loadtxt(PATH_BASE.joinpath("Results", "OpenOpen", "Test02", "inputs", "prot_xy_init.dat")).reshape(7,13,2)
    lo = np.array([-0.2]*7)
#   hi = np.array([0.1, 0.2, 0.4, 0.5, 0.6, 0.65])
#   hi = np.array([0.4, 1.75, 2.50, 3.05, 3.55, 3.95])
    hi = np.array([6, 5.55, 5.75, 5.95, 6.25, 7.85])

    fig, ax = plt.subplots(2,3)
    ax = ax.reshape(ax.size)
    for i in range(6):
        j = i
        k = i + 1
#       idx = ((etot[:,j] - etot[:,k]) < lo[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
#       plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, idx, ec=ec[j])
#       print(i, np.sum(idx))

        idx = ((etot[:,j] - etot[:,k]) > hi[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[i], clig[k], cprot[:,k], seq, idx, ec=ec[k])
        print(i, np.sum(idx))


def no_def_energy(protseq, ligseq, angle, energy):
    pair_energy = np.array([energy[int(protseq[i])-1, int(ligseq[i])-1] for i in [0,1,2]])
    if abs(angle - np.pi / 3) < 0.001:
        return - np.sum(pair_energy)
    elif angle > np.pi / 3:
        chord = 2 * np.sin((np.pi / 3 - angle)/2)
        edx =  np.exp(-chord**2/0.3**2)
        return -(pair_energy[1] + max(pair_energy[[0,2]]) + min(pair_energy[[0,2]]) * edx)
    else:
        return -pair_energy[[0,2]].max()


def get_no_def_energy(seq, angles, emat):
    return np.array([[AT.no_def_energy(s[-3:], '222', a, emat) for s in seq] for a in angles]).T


def plot_gap_03(etot, angles, enodef):
    fig, ax = plt.subplots()
    for i in range(int(etot.shape[1]/2)):
        lbl = '_'.join([str(round(x,1)) for x in angles[i*2:(i+1)*2]])
        idx = (etot[:,i*2] < enodef[:,i*2]) & (etot[:,i*2+1] < enodef[:,i*2+1])
        gap = etot[idx,i*2] - etot[idx,i*2+1]
        sns.distplot(gap, label=lbl, ax=ax)

    ax.legend()
    ax.set_xlabel("Recognition gap")


def plot_all_best_03(clig, cprot, seq, etot, enodef):
    ec = np.loadtxt(PATH_BASE.joinpath("Results", "SymAsym", "Test01", "inputs", "prot_xy_init.dat")).reshape(12,13,2)
    op = [np.less]*4 + [np.greater]
    lim = [-1.9, -1.5, -0.95, -0.15, 0.15]

    fig, ax = plt.subplots(2,5)
    for i, ii in enumerate([0,1,2,3,5]):
        j = ii*2
        k = ii*2 + 1
        idx = (op[i]((etot[:,j] - etot[:,k]), lim[i])) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, idx, ec=ec[j])
        print(i, np.sum(idx))

        plot_best_01(fig, ax[1,i], clig[k], cprot[:,k], seq, idx, ec=ec[k])


def plot_all_best_04(clig, cprot, seq, etot, enodef):
    ec = np.loadtxt(PATH_BASE.joinpath("Results", "ClampOpenSymAsym", "Test02", "inputs", "prot_xy_init.dat")).reshape(10,13,2)
#   lo = np.array([-0.01, -0.005, -0.015, -0.02, -0.02])
#   hi = np.array([0.02, 0.06, 0.08, 0.1, 0.13])
    lo = np.array([-0.08, -0.050, -0.15, -0.24, -0.22])
    hi = np.array([0.08, 0.31, 0.35, 0.46, 0.68])

    fig, ax = plt.subplots(2,5)
    for i in range(5):
        j = i*2
        k = i*2 + 1
        idx = ((etot[:,j] - etot[:,k]) < lo[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, idx, ec=ec[j])
        print(i, np.sum(idx))

        idx = ((etot[:,j] - etot[:,k]) > hi[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[1,i], clig[k], cprot[:,k], seq, idx, ec=ec[k])
        print(i, np.sum(idx))

        X = (etot[:,j] - etot[:,k])[(enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])]
        print(X.min(), X.max())


def check_triangle_orientation(ec, cprot, adj):
    triangles = utils.find_triangles(adj)
    orient0 = utils.wedge_product(ec[triangles])
    orient1 = utils.wedge_product(cprot[:,:,triangles])
    return orient1.T.dot(orient0) == len(triangles)


#########################################################
        ### Identify recognition


def update_energy_and_save(path, etot, ent, is_valid, enodef):
    etot = np.min([etot, enodef], axis=0)
    etot[is_valid==False] = np.nan
    ftot = (etot.T + ent).T
    np.save(path.joinpath("energy.npy"), etot)
    np.save(path.joinpath("free_energy.npy"), ftot)


### For the MeCh model we looked at 34 ligand shapes, in which case "t=1"
### For the G-MeCh model we looked at 14 ligand shapes, in which case "t=0"
def process_test_01(path, t=0):
    print(path)
    clig, cprot, etot, edef, echem, einit, lseq, angles, kmat, emat, is_valid, enodef, ent, seq = load_test_results(path)
    update_energy_and_save(path, etot, ent, is_valid, enodef)
    ftot = (etot.T + ent).T

    cols = ["ID", "no", "pseq", "kmat", "emat", "ent", "detot", "dedef", "dechem"] + \
           [f"{a}{b}" for a in ["lseq", "ang", "clig", "cprot", "etot", "edef", "echem", "enodef", "ftot"] for b in [1,2]]
           
    df0 = pd.DataFrame(columns=cols)
    df_list = []
    for i in range([14,34][t]):
#   for i in range(14):
        if t == 0:
            i1 = i * 2
            i2 = i * 2 + 1
        elif t == 1:
            i1 = i
            i2 = i + 1

        # Label and count valid / invalid outputs
        idx = (etot[:,i1] < enodef[:,i1]) & (etot[:,i2] < enodef[:,i2]) & is_valid[:,i1] & is_valid[:,i2]
        N = np.sum(idx)
        print(f"# invalid configs: {np.sum((is_valid[:,i1]==False) | (is_valid[:,i2]==False))}")
        print(f"# low energy: {np.sum((etot[:,i1] >= enodef[:,i1]) | (etot[:,i2] >= enodef[:,i2]))}")
        if not N:
            continue

        df = df0.copy()
        df['ID'] = [path.parts[-2]] * N
        df['no'] = [path.parts[-1]] * N
        df['pseq'] = list(seq[idx])
        df['ent'] = ent[idx]
        df['kmat'] = [kmat] * N
        df['emat'] = [emat] * N
        df['lseq1'] = [lseq[i1]] * N
        df['lseq2'] = [lseq[i2]] * N
        df['ang1'] = [angles[i1]] * N
        df['ang2'] = [angles[i2]] * N
        df['clig1'] = [clig[i1]] * N
        df['clig2'] = [clig[i2]] * N
        for k, d in zip(["cprot", "etot", "edef", "echem", "enodef", "ftot"], [cprot, etot, edef, echem, enodef, ftot]):
            for j in [0,1]:
                df[f"{k}{j+1}"] = list(d[idx,[i1,i2][j]])
            if k in ["etot", "edef", "echem", "enodef"]:
                df[f"d{k}"] = d[idx,i2]  - d[idx,i1]
        df_list.append(df)
    df = pd.concat(df_list, ignore_index=True)

    adj = utils.break_links(utils.get_adj_from_coord(load_equil_coord(path)))
    df['gap'] = np.abs(df['detot'])
    df['dE'] = utils.bound_dechem_df(df)
    df['meanE'] = df.pseq.apply(lambda x: np.sum([emat[int(s)-1,1] for s in x[-3:]]))
    df['meanK'] = df.pseq.apply(lambda x: np.sum(GE.getK(adj, kmat, x[:-3])) / np.sum(adj) )
    df['krat'] = df.kmat.apply(lambda x: x.max() / x.min())
    df['erat'] = df.emat.apply(lambda x: x.max() / x.min())
    return df


def make_histograms(df):
    bins, combos = utils.get_plot_bins()
    out = {}
    for c1, c2 in combos:
        xbin, ybin = bins[c1], bins[c2]
        hist = np.histogram2d(df[c1].values, df[c2].values, bins=[xbin, ybin])[0]
        out[f"{c1}_{c2}"] = {'bins':(xbin, ybin), 'hist':hist}
    return out


def add_all_histograms(base):
    bins, combos = utils.get_plot_bins()
    lbls = ["all", "bind1", "bind2", "rec1", "rec2"]
    out = {}
    for l in lbls:
        files = list(base.glob(f"*/hist_{l}.pickle"))
        data = pickle.load(open(files[0], 'rb'))
        
        out[l] = data.copy()
        for f in files[1:]:
            data = pickle.load(open(f, 'rb'))
            for c1, c2 in combos:
                out[l][f"{c1}_{c2}"]['hist'] = out[l][f"{c1}_{c2}"]['hist'] + data[f"{c1}_{c2}"]['hist']
    pickle.dump(out, open(base.joinpath("hist_summary.pickle"), 'wb'))
    


def summarize_one(path, cut=2):
    df = pd.read_pickle(path)
    cols = ["detot", "dedef", "dechem", "meanE", "meanK", "krat", "erat"] + \
           [f"{a}{b}" for a in ["ang", "etot", "edef", "echem", "ftot"] for b in [1,2]]
    mean = df.loc[:, cols].mean()
    std = df.loc[:, cols].std()

    i1 = (df.ftot1<0)&(df.detot>cut)
    i2 = (df.ftot2<0)&(df.detot<-cut)

    i3 = (df.ftot1<0)&(df.ftot2>-cut)&(df.detot>cut)
    i4 = (df.ftot2<0)&(df.ftot1>-cut)&(df.detot<-cut)

    idx_list = [i1, i2, i3, i4]
    print(path)
    print([np.sum(i) for i in idx_list])

    lbls = ["all", "bind1", "bind2", "rec1", "rec2"]
    df_list = [df.loc[i] for i in idx_list]
    for df0, l in zip([df] + df_list, lbls):
        out = make_histograms(df0)
        pickle.dump(out, open(Path(path).with_name(f"hist_{l}.pickle"), 'wb'))

    return path, [cols, mean, std], idx_list


def summarize_all():
    files = sorted(glob.glob("*/processed.pickle"))
    with Pool(N_PROC) as pool:
        res = pool.map(summarize_one, files, 11)
    cols = res[0][1][0]
    mean = np.array([r[1][1] for r in res])
    std = np.array([r[1][2] for r in res])
    pickle.dump([cols, mean, std], open("stats.pickle", "wb"))

    lbls = ["bind1", "bind2", "rec1", "rec2"]
#   data = [[], [], [], []]
#   for r in res:
#       df0 = pd.read_pickle(r[0])
#       for i, l in enumerate(lbls):
#           data[i].append(df0.loc[r[2][i]])

    for i, l in enumerate(lbls):
        recog = pd.concat([pd.read_pickle(r[0]).loc[r[2][i]] for r in res], ignore_index=True, axis=0)
#       recog = pd.concat(data[i], ignore_index=True, axis=0)
        print(l, len(recog))
        try:
            recog.to_pickle(f"{l}.pickle")
        except Exception as e:
            print(e)

    add_all_histograms(Path.cwd())


def make_figures():
    df_rec = pd.read_pickle("recognition.pickle")
#   df_all = pd.read_pickle("all_res.pickle")
#   df_list = [df_all, df_rec, df_rec.loc[df_rec.edef2>0.1]]
    df_list = [df_rec, df_rec.loc[df_rec.edef2>0.1]]
#   lbls = ["all", "rec-all", "rec-def"]
    lbls = ["rec-all", "rec-def"]
    cols = ["detot", "dedef", "dechem", "meanE", "meanK", "krat", "erat"] + \
           [f"{a}{b}" for a in ["ang", "etot", "edef", "echem", "ftot"] for b in [1,2]]

    combos = [('ent', 'detot'),
              ('ent', 'dechem'),
              ('ent', 'dedef'),
              ('meanK', 'detot'),
              ('meanE', 'detot'),
              ('dedef', 'dechem'),
              ('dedef', 'detot'),
              ('dechem', 'detot'),
              ('echem2', 'entdef2')]

    base = Path("/data/Jmcbride/RecProt/ClampOpen/Figures/")

    for df, l in zip(df_list, lbls):
        df['entdef2'] = df['ent'] + df['edef2']
        for c in cols:
            fig, ax = plt.subplots()
            sns.distplot(df[c], ax=ax)
            fig.savefig(base.joinpath(f"1d_{c}_{l}.png"), bbox_inches='tight')
            plt.close()

        for c1, c2 in combos:
            fig, ax = plt.subplots()
            hist, xedges, yedges = np.histogram2d(df[c1].values, df[c2].values, bins=20)
            plots.plot_heatmap(hist, xedges, yedges, c1, c2, fig, ax)
            fig.savefig(base.joinpath(f"2d_{c1}_{c2}_{l}.png"), bbox_inches='tight')
            plt.close()
            

def run_tests_fig3():
    angles = np.arange(2.5, 90, 2.5) * np.pi / 180
    degrees = 2 * (90 - angles * 180 / np.pi)
    ener = np.load("free_energy.npy")

    i0 = [17]
    di_arr = np.arange(-1, -9, -1)
    
        



    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", action="store", default="process", type=str)
    parser.add_argument("-i", action="store", default=0, type=int)
    parser.add_argument("-t", action="store", default=0, type=int)
    return parser.parse_args()





if __name__ == "__main__":
#   angle_pairs_equal_deformation()
    
#   write_test_ligands_01()
#   write_test_ligands_02()

    args = parse_args()

    if args.mode == "process":
        path = Path.cwd()
        df = process_test_01(path, t=args.t)
        df.to_pickle("processed.pickle")


    elif args.mode == "compile":
        summarize_all()

    elif args.mode == "plot":
        make_figures()


