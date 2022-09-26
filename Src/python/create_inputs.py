"""
Tools for creating input files for main FORTRAN program
"""

from itertools import permutations, product
import os
from pathlib import Path
import pickle
import shutil
import time

import numpy as np

from rec_io import PATH_BASE
import rec_io
import rec_utils as utils



# Define a basic lattice for creating protein shapes
def create_hex_lattice(nx, ny):
    nPart = nx*ny*2
    xy = np.zeros( (nPart,2), dtype=float)

    dy = 0.5 * 3**0.5

    for j in range(ny):
        for i in range(nx):
            xy[j*(nx*2)+i*2,0] = i
            xy[j*(nx*2)+i*2,1] = (j*2) * dy

            xy[j*(nx*2)+i*2+1,0] = i + 0.5
            xy[j*(nx*2)+i*2+1,1] = (j*2+1) * dy

    # Remove one column to create x-symmetry
    xy = xy[(xy[:,0] > xy[:,0].min())]
    
    # Reorder indices so that they can be read from left to right, bottom to top
    idx = np.argsort(xy[:,0] + xy[:,1] * 20)
    xy = xy[idx]

    xy[:,0] = xy[:,0] - xy[:,0].mean()
    xy[:,1] = xy[:,1] - xy[:,1].min() - 0.5 * 3**0.5

    return xy


# Define a basic set of protein shapes to work with
def protein_shapes():
    basic_lattice = create_hex_lattice(9, 6)
    shapes = [{'N':7, 'name':'w1h1', 'pidx':[3, 4, 11, 12, 13, 20, 21], 'bidx':[0,3,1]},
              {'N':13, 'name':'w2h1', 'pidx':[2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 20, 21, 22],
               'bidx':[1, 6, 2]},
              {'N':19, 'name':'w3h1',
               'pidx':[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23],
               'bidx':[2, 9, 3]},
              {'N':25, 'name':'w4h1',
               'pidx':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
               'bidx':[3, 12, 4]},
              {'N':16, 'name':'w2h2-',
               'pidx':[2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 20, 21, 22, 28, 29, 30],
               'bidx':[1, 6, 2]},
              {'N':18, 'name':'w2h2',
               'pidx':[2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 20, 21, 22, 27, 28, 29, 30, 31],
               'bidx':[1, 6, 2]},
              {'N':18, 'name':'w2h3-',
               'pidx':[2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 20, 21, 22, 28, 29, 30, 37, 38],
               'bidx':[1, 6, 2]},
              {'N':22, 'name':'w2h3',
               'pidx':[2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 20, 21, 22, 27, 28, 29, 30, 31, 36, 37, 38, 39],
               'bidx':[1, 6, 2]},
              {'N':24, 'name':'w3h2-',
               'pidx':[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 27, 28, 29, 30, 31],
               'bidx':[2, 9, 3]},
              {'N':26, 'name':'w3h2',
               'pidx':[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 30, 31, 32],
               'bidx':[2, 9, 3]},
              {'N':28, 'name':'w3h3--',
               'pidx':[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 27, 28, 29, 30, 31, 36, 37, 38, 39],
               'bidx':[2, 9, 3]},
              {'N':32, 'name':'w3h3-',
               'pidx':[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 30, 31, 32, 35, 36, 37, 38, 39, 40],
               'bidx':[2, 9, 3]}
             ]

    for i in range(len(shapes)):
        xy = basic_lattice[shapes[i]['pidx']]
        links = [[shapes[i]['bidx'][j] for j in [0,2]]]
        adj = utils.break_links(utils.get_adj_from_coord(xy), links)
        shapes[i].update({'xy':xy, 'adj':adj})

    return shapes


# Get a sequence list for a given sequence length.
# If the total number of sequences exceeds "cut", sample without replacement.
# In the final work, we did not report any sequences longer than "cut"
def pick_seq(l, cut=2**16):
    if l <= 25:
        seq = np.loadtxt(PATH_BASE.joinpath("Inputs", "protein_seq", f"amino{l}.dat"))
        if cut == 0:
            return seq
        elif len(seq) > cut:
            seq = [''.join(s.astype(int).astype(str)) for s in seq]
            return np.array([list(s) for s in np.random.choice(seq, size=cut, replace=False)])
        else:
            return seq
    else:
        return generate_protein_seq(l, cut)
    

def write_protein_shapes():
    shapes = protein_shapes()
    for s in shapes:
        if s['N'] > 22:
            continue
        print(s['name'], s['N'])
        base = PATH_BASE.joinpath("Inputs", "protein_shapes", s['name'])
        base.mkdir(parents=True, exist_ok=True)
        np.savetxt(base.joinpath("prot_xy.dat"), s['xy'])

        adj = np.array([[i,j] for i, j in zip(*np.where(s['adj'])) if j > i]) + 1
        np.savetxt(base.joinpath("prot_adjacency.dat"), adj, fmt='%d')

        seq = pick_seq(s['N']+len(s['bidx']), 0).astype(int)
        np.savetxt(base.joinpath("prot_seq.dat"), seq, fmt='%d')

        params = rec_io.init_params()
        new_vals = {'nprot':len(seq), 'namino':seq.shape[1]-3, 'plinks': len(adj),
                    'ibind':[i+1 for i in s['bidx']]}
        params = rec_io.update_param_vals(params, new_vals)
        rec_io.write_params(base.joinpath("parameters.txt"), params)


def get_protein_seq(l):
    return product(['1', '2'], repeat=l)


# A quick way of sampling without replacement is to generate extra samples,
# and only use the first N unique samples
def generate_protein_seq(l, N, replace=False):
    if replace:
        return np.random.randint(2, size=N*l).reshape(N,l).astype(str)
    else:
        seq = set([''.join((x+1).astype(str)) for x in np.random.randint(2, size=N*l*2).reshape(2*N,l)])
        return np.array([list(s) for s in list(seq)[:N]])


def write_protein_sequences(lo=6, hi=23):
    for i in range(lo, hi):
        with open(PATH_BASE.joinpath("Inputs", "protein_seq", f"amino{i}.dat"), 'w') as o:
            ts = time.time()
            for a in get_protein_seq(i):
                o.write(f"{' '.join(a)}\n")
            print(f"{(time.time()-ts)/60}")



def basic_ligand(ang, r=1):
    x = [-np.cos(ang)*r, 0, np.cos(ang)*r]
    y = [-np.sin(ang)*r, 0, -np.sin(ang)*r]
    return np.array([x, y]).T


def update_prot_ypos(cprot, i0, ligands, extra=0.001):
    out = []
    for i, l in enumerate(ligands):
        clig = l[2]
        if cprot[i0,0] <= clig[0,0]:
            ypos = clig[0,1] + extra
        else:
            dx = (clig[1,1] - clig[0,1])
            dy = (clig[1,0] - clig[0,0])
            if dx == 0:
                y = clig[0,1] + extra
            else:
                dxdy = dx / dy
                y = clig[0,1] + dxdy * (cprot[i0,0] - clig[0,0]) + extra

        c = cprot.copy()
        c[:,1] = c[:,1] - c[0,1] + y
        out.append(c.flatten())
    return np.array(out)


def create_big_small_ligand(angle=np.pi/3):
    R = np.arange(0.4, 1.7, 0.1)
    ligand_colors = [2, 2, 2]
    ligands = []
    for r in R:
        ligands.append((ligand_colors, angle, basic_ligand(angle, r)))
    rec_io.write_ligands(PATH_BASE.joinpath("Inputs", "ligand_sets", "bigsmall13.dat"), ligands)


### We didn't use these ligands (R!=1) in the final work,
### as this type of binding required extensive stretching of single bonds,
### and collective motion cannot help binding. 
def create_all_ligand():
    R = np.arange(0.4, 1.7, 0.1)
    angles = np.arange(30, 95, 5) * np.pi / 180
    ligand_colors = [2, 2, 2]
    ligands = []
    for r in R:
        for a in angles:
            ligands.append((ligand_colors, a, basic_ligand(a, r)))
    rec_io.write_ligands(PATH_BASE.joinpath("Inputs", "ligand_sets", "allligand13.dat"), ligands)




if __name__ == "__main__":

    write_protein_sequences()
    write_protein_shapes()
    create_big_small_ligand()
    create_all_ligand()


 
