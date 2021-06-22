from itertools import permutations, product
import os
from pathlib import Path
import pickle

import numpy as np

from rec_io import PATH_BASE



def get_protein_seq(i):
    return product(['1', '2'], repeat=i)


def write_protein_sequences():
    for i in range(6, 19):
        with open(PATH_BASE.joinpath("Inputs", f"amino{i}.dat"), 'w') as o:
            for a in get_protein_seq(i):
                o.write(f"{' '.join(a)}\n")


def write_ligand_sequences(path, ligands):
    with open(path, 'w') as o:
        for l in ligands:
            o.write("{0:d} {1:d} {2:d} {3:f}\n".format(*l))



if __name__ == "__main__":

    write_protein_sequences()

 
