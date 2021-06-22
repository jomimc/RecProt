from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import numpy as np
from palettable.scientific.diverging import Berlin_20, Roma_20
from scipy.signal import argrelmin
import seaborn as sns

import create_inputs
import rec_io
from rec_io import PATH_BASE



#########################################################
### Load initial results

### These came from running 180 ligands that only differ
### by angle. Only one protein spring, and only one
### chemical bond strength.


def load_results(path):
    data = rec_io.load_output(path)
    etot, edef, echem, einit = data[:,9:13].T
    clig = data[:,2:8]
    cprot = data[:,13:]
    seq = data[:,8]
    return clig, cprot, etot, edef, echem, einit, seq


def load_equil_coord():
    return np.loadtxt(PATH_BASE.joinpath("Src", "test_angle", "inputs", "prot_xy.dat"))


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





#########################################################
### Analyse test results


def load_test_results(p1="ClampOpen", p2="Test02", nlig=18):
    path = PATH_BASE.joinpath("Results", p1, p2, "output.dat")
    clig, cprot, etot, edef, echem, einit, seq = load_results(path)
    clig = clig[:nlig]
    cprot = cprot.reshape(2**16, nlig, 26)
    etot = etot.reshape(2**16,nlig)
    edef = edef.reshape(2**16,nlig)
    echem = echem.reshape(2**16,nlig)
    einit = einit.reshape(2**16,nlig)
    seq = np.array([list(str(int(s))) for s in seq[::nlig]], float)

    angles = np.loadtxt(PATH_BASE.joinpath("Results", p1, p2, "inputs", "ligands.dat"))[:,3]
    return clig, cprot, etot, edef, echem, einit, angles, seq


def plot_gap_01(etot, angles, enodef):
    fig, ax = plt.subplots(2,1)
    for i in range(int(etot.shape[1]/2)):
        lbl = '_'.join([str(round(x,1)) for x in angles[i*2:(i+1)*2]])
        idx = (etot[:,i] < enodef[:,i*2]) & (etot[:,i+1] < enodef[:,i*2+1])
        gap = etot[idx,i*2] - etot[idx,i*2+1]
        j = i//5
        sns.distplot(gap, label=lbl, ax=ax[j])

    for a in ax:
        a.legend()
        a.set_xlabel("Recognition gap")


def plot_best_01(fig, ax, clig, cprot, seq, idx, ec=''):
    ax.plot(*ec.T, 'sk', fillstyle='none', ms=10)
    if isinstance(ec, str) or 1:
        ec = load_equil_coord()
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
        lo = np.array([-0.015, -0.02, -0.05, -0.1])
        hi = np.array([0.2, 1, 2, 2])
#   lo, hi = lo * 10, hi * 10

    fig, ax = plt.subplots(2,itot)
    for i in range(itot):
        j = istart + i*2
        k = istart + i*2 + 1
        plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, (etot[:,j] - etot[:,k]) < lo[i], ec=ec[j])
        plot_best_01(fig, ax[1,i], clig[k], cprot[:,k], seq, (etot[:,j] - etot[:,k]) > hi[i], ec=ec[k])


def plot_gap_02(etot, angles, enodef):
    fig, ax = plt.subplots()
    for i in range(etot.shape[1]-1):
        lbl = '_'.join([str(round(x,2)) for x in angles[i:i+2]])
        idx = (etot[:,i] < enodef[:,i]) & (etot[:,i+1] < enodef[:,i+1])
        gap = etot[idx,i] - etot[idx,i+1]
        sns.distplot(gap, label=lbl)
    ax.legend()
    ax.set_xlabel("Recognition gap")


def plot_all_best_02(clig, cprot, seq, etot, enodef):
    ec = np.loadtxt(PATH_BASE.joinpath("Results", "OpenOpen", "Test02", "inputs", "prot_xy_init.dat")).reshape(7,13,2)
    lo = np.array([-0.2]*7)
#   hi = np.array([0.1, 0.2, 0.4, 0.5, 0.6, 0.65])
    hi = np.array([0.4, 0.9, 2.3, 2.7, 3.0, 3.5])

    fig, ax = plt.subplots(2,6)
    for i in range(6):
        j = i
        k = i + 1
        idx = ((etot[:,j] - etot[:,k]) < lo[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[0,i], clig[j], cprot[:,j], seq, idx, ec=ec[j])
        print(i, np.sum(idx))

        idx = ((etot[:,j] - etot[:,k]) > hi[i]) & (enodef[:,j] > etot[:,j]) & (enodef[:,k] > etot[:,k])
        plot_best_01(fig, ax[1,i], clig[k], cprot[:,k], seq, idx, ec=ec[k])
        print(i, np.sum(idx))


def no_def_energy(protseq, ligseq, angle, energy):
    pair_energy = np.array([energy[int(protseq[i])-1, int(ligseq[i])-1] for i in [0,1,2]])
    if angle > np.pi / 3:
        return -max(pair_energy[:2].sum(), pair_energy[1:].sum())
    else:
        return -pair_energy[[0,2]].max()


def get_no_def_energy(seq, angles):
    emat = np.array([[1, 1.5], [1.5, 2]])
    return np.array([[AT.no_def_energy(s[13:], '222', a, emat) for s in seq] for a in angles]).T



if __name__ == "__main__":
#   angle_pairs_equal_deformation()
    
    write_test_ligands_01()
    write_test_ligands_02()

    pass


