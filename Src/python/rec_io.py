import os
from pathlib import Path
import pickle

import numpy as np

PATH_BASE = [p for p in Path.cwd().parents if p.stem == 'RecProt'][0]


def load_output(path):
    return np.loadtxt(path, skiprows=2)



