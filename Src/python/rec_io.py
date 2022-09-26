import os
from pathlib import Path
import pickle

import numpy as np

import rec_utils as utils

PATH_BASE = [p for p in Path.cwd().parents if p.stem == 'RecProt'][0]

PARAM_KEYS = ["nprot", "namino", "plinks", "nligand", "nbind", "ibind",
              "ks", "kw", "kww", "kstep", "kmod", "els", "elw", "lstruct",
              "dr11", "dr12", "dr22", "ltranslation"]


def load_output(path):
    return np.loadtxt(path, skiprows=2)


def parse_params(path):
    params = {}
    name = ''
    for l in open(path, 'r'):
        if len(name):
            if name == 'ibind':
                params[name] = [int(x) for x in l.strip('\n').split()]
            else:
                params[name] = utils.eval_input(l.strip('\n'))
            name = ''
        else:
            name = l.strip('\n').split('::')[1].replace(' ', '')
    return params


def get_spring_const(params):
    return np.array([[params['kww'], params['kw']], 
                     [params['kw'], params['ks']]])


def get_chem_const(params):
    ks, kw = params['els'], params['elw']
    km = 0.5 * (kw + ks)
    return np.array([[kw, km], [km, ks]])


def parse_ligands(path):
    data = np.loadtxt(path)
    N = int((data.shape[1] - 1)/3)
    seq = data[:,:N].astype(int).astype(str)
    angles = data[:,N].astype(float)
    x = data[:,N+1:N*2+1]
    y = data[:,N*2+1:N*3+1]
    xy = np.array([x.flatten(), y.flatten()]).T.reshape(len(data),N,2)
    return seq, angles, xy


def init_params():
    text = ["# Number of proteins :: nprot  :: Number of protein sequences",
            "# Protein size :: namino  :: Sequence length (N positions, and 3 binding sites)",
            "# Protein bonds :: plinks  :: Number of internal protein bonds",
            "# Number of ligands ::  nligand",
            "# Ligand size ::  nbind",
            "# Protein binding sites :: ibind :: Indices of protein binding sites",
            "# Elastic energy :: ks  ::  energy constant for stiff springs",
            "# Elastic energy :: kw  ::  energy constant for medium springs",
            "# Elastic energy :: kww  ::  energy constant for weak springs",
            "# Elastic energy :: kstep  ::  number of iterations using different spring constants",
            "# Elastic energy :: kmod  :: scalar modifier for ramping up spring constant",
            "# Electrostatic energy :: els  :: energy constant for strong chemical interaction",
            "# Electrostatic energy :: elw  :: energy constant for weak chemical interaction",
            "# Use local structure? :: lstruct  :: If no, all sequences have same structure",
            "# Structure variation :: dr11  :: magnitude / direction of variation",
            "# Structure variation :: dr12  :: magnitude / direction of variation",
            "# Structure variation :: dr22  :: magnitude / direction of variation",
            "# Translation step :: lTranslation :: the protein first translates downward before deforming"]

    default_vals = [65536, 13, 25, 10, 3, [2, 7, 3], 20, 8, 2, 1, 1, 2, 1, "F", 0., 0., 0., "T"]

    return {k:{'txt':t, 'val':v} for k, t, v in zip(PARAM_KEYS, text, default_vals)}


def update_param_vals(old, new):
    for k, v in new.items():
        old[k]['val'] = v
    return old


def write_params(path, params):
    with open(path, 'w') as o:
        for k in PARAM_KEYS:
            o.write(f"{params[k]['txt']}\n")
            if isinstance(params[k]['val'], (float, str, int)):
                o.write(f"{params[k]['val']}\n")
            elif isinstance(params[k]['val'], (list, np.ndarray)):
                s = ' '.join([f"{v}" for v in params[k]['val']])
                o.write(f"{s}\n")
            else:
                print(f"Wrong type for {k}: params[k]['val']")


def write_ligands(path, ligands):
    with open(path, 'w') as o:
        for seq, ang, xy in ligands:
            s_seq = ' '.join([f"{s:d}" for s in seq])
            s_ang = f"{ang:f}"
            s_xy  = ' '.join([f"{x:f}" for x in xy.T.flatten()])
            o.write(f"{s_seq} {s_ang} {s_xy}\n")




