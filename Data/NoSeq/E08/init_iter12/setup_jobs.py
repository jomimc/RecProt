from copy import deepcopy
import os
from pathlib import Path
import pickle
import shutil
import sys

sys.path.insert(0, "/home/jmcbride/RecProt/Src/python")

import numpy as np

import rec_io


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        p = Path(d)
        p.mkdir(exist_ok=True, parents=True)

        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def const_set():
    Kmax = np.round(10.**np.linspace(1, 5, 300), 2)

    for ks in Kmax:
        yield {"ks":ks, "kw": ks, "kww": ks}


#ef write_params(path, params):
#   lbls = {'Elastic energy':['ks', 'kw', 'kww'], 'Electrostatic energy':['els', 'elw']}
#   lbls_list = [f"# {a} :: {b}" for a, l in lbls.items() for b in l]
#   params_list = [f"{p:09.3f}" for p in params]
#   with open(path, 'w') as o:
#       for l, p in zip(lbls_list, params_list):
#           o.write(f"{l}\n")
#           o.write(f"{p}\n")
        

def setup_one(base, params):
    template = "Template"
#   if not os.path.exists(base):
#       copytree(template, base)
    path = os.path.join(base, "inputs", "parameters.txt")
    rec_io.write_params(path, params)


def setup_all():
    default_params = rec_io.init_params()
    update_params = rec_io.parse_params("Template/inputs/parameters.txt")
    params = rec_io.update_param_vals(default_params, update_params)
    param_list = []
    for i, const in enumerate(const_set()):
        params = rec_io.update_param_vals(params.copy(), const)
        param_list.append(deepcopy(params))
        base = f"{i:03d}"
#       shutil.rmtree(base)
        setup_one(base, params)
    pickle.dump(param_list, open("param_list.pickle", "wb"))


if __name__ == "__main__":

    setup_all()
        

