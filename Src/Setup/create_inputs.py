from itertools import permutations, product
import os
from pathlib import Path
import pickle

import numpy as np


PATH_BASE = [p for p in Path.cwd().parents if p.stem == 'RecProt'][0]


def get_protein_seq(i):
    return product(['1', '2'], repeat=i)


def write_protein_sequences():
    for i in range(6, 19):
        with open(PATH_BASE.joinpath("Inputs", f"amino{i}.dat"), 'w') as o:
            for a in get_protein_seq(i):
                o.write(f"{' '.join(a)}\n")


if __name__ == "__main__":

    write_protein_sequences()

 
