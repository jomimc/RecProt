from pathlib import Path

import numpy as np
from scipy.signal import argrelmin

import create_inputs
import rec_io
from rec_io import PATH_BASE


def load_results():
    data = rec_io.load_output(PATH_BASE.joinpath("Results", "test_angle", "output.dat"))
    etot, edef, echem = data[:,10:13].T
    cprot = data[:,13:]
    return cprot, etot, edef, echem


def angle_pairs_equal_deformation():
    cprot, etot, edef, echem = load_results()
    deform = np.arange(0.1, 0.6, 0.1)
    total = [-4.4, -4.2, -4.0, -3.8]
    angles = np.arange(180) * np.pi / 2 / 180

#   for d in deform:
#       idx = argrelmin(np.abs(edef-d))
#       print(f"For e_def = {d}, closest match is found for...")
    for d in total:
        idx = argrelmin(np.abs(etot-d))
        print(f"For e_tot = {d}, closest match is found for...")
        print("Deformation energy: ", edef[idx])
        print("Total energy: ", etot[idx])
        print("Angles: ", (angles*180/np.pi)[idx])
        print("")


def test_ligands_to_write():
    cprot, etot, edef, echem = load_results()
    deform = np.arange(0.1, 0.6, 0.1)
    total = [-4.4, -4.2, -4.0, -3.8]
    angles = np.arange(180) * np.pi / 2 / 180

    ligands = []
    ligand_colors = [2, 2, 2]

    for d in deform:
        idx = argrelmin(np.abs(edef-d))
        for a in angles[idx]:
            ligands.append(ligand_colors + [a])

    for d in total:
        idx = argrelmin(np.abs(etot-d))
        for a in angles[idx]:
            ligands.append(ligand_colors + [a])
    return ligands


def write_test_ligands():
    ligands = test_ligands_to_write()
    path = PATH_BASE.joinpath("Results", "Template", "ligands.dat")
    create_inputs.write_ligand_sequences(path, ligands)


if __name__ == "__main__":
#   angle_pairs_equal_deformation()
    
    write_test_ligands()



