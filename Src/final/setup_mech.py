from pathlib import Path
import shutil

import numpy as np
from scipy.signal import argrelmin

import angle_test as AT
import create_inputs as CI
import rec_io
from rec_io import PATH_BASE
import rec_utils as utils


def all_angles():
    clig, cprot, etot, edef, echem, einit, seq = AT.load_results(PATH_BASE.joinpath("Src", "test_angle", "output.dat"))
    new_angles = np.arange(0, 185, 5) * np.pi / 2 / 180
    angles = np.arange(180) * np.pi / 2 / 180

    ligands = []
    init_config = []
    ligand_colors = [2, 2, 2]

    # Controlling for shape mismatch
    for a in new_angles:
        idx = argrelmin(np.abs(angles-a))
        for i in idx[0]:
            ligands.append((ligand_colors, angles[i], CI.basic_ligand(angles[i])))
            c = cprot[i]
            c[1::2] = c[1::2] + 0.001
            init_config.append(c)
        
    return ligands, np.array(init_config)


def write_test_ligands_01(update=0):
    ligands, cprot = all_angles()
    shape = CI.protein_shapes()[1]

    base = PATH_BASE.joinpath("Inputs", "protein_shapes", shape['name'])
    dest = PATH_BASE.joinpath("Results", "NoSeq",  shape['name'], "Template")
    dest.joinpath("inputs").mkdir(parents=True, exist_ok=True)

    shutil.copy(PATH_BASE.joinpath("Src", "fortran", "protein.exe"), dest.joinpath("protein.exe"))

    params = rec_io.init_params()
    params = rec_io.update_param_vals(params, {'nprot':1, 'nligand':len(ligands), 'els':8, 'elw':8})
    rec_io.write_params(dest.joinpath("inputs", "parameters.txt"), params)

    seq = np.ones((1,shape['N'] + len(shape['bidx'])), int) * 2
    np.savetxt(dest.joinpath("inputs", "prot_seq.dat"), seq, fmt='%d')

    rec_io.write_ligands(dest.joinpath("inputs", "ligands.dat"), ligands)
    shutil.copy(base.joinpath("prot_xy.dat"), dest.joinpath("inputs", "prot_xy.dat"))
    shutil.copy(base.joinpath("prot_adjacency.dat"), dest.joinpath("inputs", "prot_adjacency.dat"))

#   np.savetxt(dest.joinpath("inputs", "prot_xy_init.dat"), cprot)
    np.savetxt(dest.joinpath("inputs", "prot_xy_init.dat"), np.array([shape['xy'].flatten()]*len(ligands)))




