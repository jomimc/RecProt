"""
Tools for studying evolvability of proteins
"""
from itertools import product
from pathlib import Path
import pickle
import time
import sys

from multiprocessing import Pool
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN


import local_utils as utils


def evolve_seq(seq_str, fit, i0=-1, N=100, T=1):
    key = {s:i for i, s in enumerate(seq_str)}
    seq_out = []
    cost = []
    if i0 < 0:
        s = np.random.choice(list(seq_str))
    else:
        s = seq_str[i0]
    j = key[s]

    for i in range(N):
        sn = mutate(s)
        k = key[sn]
        df = fit[j] - fit[k]
#       print(i, j, k, ftot1[k], ftot2[k], c, cn)

        toss = np.random.rand()
        if df <= 0:
            j = k
            s = sn
        elif toss <= np.exp(-df/T):
            j = k
            s = sn

        seq_out.append(j)
        cost.append(fit[j])
    return np.array(cost)


#ef evolve_new_fun(seq_str, fmin, i0, j0, avoid, N=100, T=1, add_avoid=False, jstart=-1, falg='loose', fcut=2):
def evolve_new_fun(path, i0, j0, avoid, jstart=-1, N=100, T=1, add_avoid=False, falg='loose', fcut=2):
    fmin = np.load(path)
    p_i = path.stem.split('_')[1]
    seq_str = np.load(f"/data/Jmcbride/RecProt/Data/Sequences/seq_str_{p_i}.npy")
    fit_i = get_fitness_map(fmin, i0, avoid)
    if add_avoid:
        if not isinstance(avoid, list):
            avoid = list(avoid)
        avoid.append(i0)
    fit_j = get_fitness_map(fmin, j0, avoid, alg=falg, cut=fcut)
#   fmin = np.min([ftot, np.zeros(ftot.shape)], axis=0)

    key = {s:i for i, s in enumerate(seq_str)}
    if jstart >= 0:
        j = jstart
    else:
        j = fit_i.argmax()
    s = seq_str[j]

    fit = []
    fi = []
    fj = []
    nfit = []

    for i in range(N):
        sn = mutate(s)
        k = key[sn]
        df = fit_j[j] - fit_j[k]

        toss = np.random.rand()
        if df <= 0:
            j = k
            s = sn
        elif toss <= np.exp(-df/T):
            j = k
            s = sn

        fi.append(fmin[j,i0])
        fj.append(fmin[j,j0])
        fit.append(fit_j[j])
        nfit.append(np.sum(fmin[j]<0))
    return [np.array(f) for f in [fit, fi, fj, nfit]]


def generate_(base, i0, j0, avoid, fit_idx, ran, rep):
    if ran:
        for i in np.random.choice(fit_idx, size=rep, replace=True):
            yield base, i0, j0, avoid, i
    else:
        for i in range(rep):
            yield base, i0, j0, avoid, -1


def get_evolution_traj_stats(seq, fmin, i0, j0, T=1, N=200, rep=500, add_avoid=False, falg='loose', fcut=2, ran=False, mp=False):
    seq_str = np.array([''.join(s.astype(int).astype(str)) for s in seq])
    avoid = np.arange(1, 12, 2)
    fit_idx = np.where(get_fitness_map(fmin, i0, avoid) > 0)[0]
    base = "/data/Jmcbride/RecProt/Data/"
    res = []
    print("Number of fit sequences: ", len(fit_idx))
    if mp:
        with Pool(60) as pool:
            res = list(pool.starmap(evolve_new_fun, generate_(base, i0, j0, avoid, fit_idx, ran, rep), 3))
    else:
        for j in range(rep):
            if ran:
                istart = np.random.choice(fit_idx)
            else:
                istart = -1
            res.append(evolve_new_fun(seq_str, fmin, i0, j0, avoid, jstart=j, N=N, T=T, add_avoid=add_avoid, falg=falg, fcut=fcut))
    return np.array(res)


def generate_2_(path, i0, j0, avoid, fit_idx, N, rep):
    Tarr = 10.**np.arange(-1, 2.5, 0.5)
    alg = ['looser', 'loose', 'strict']
    add_avoid = [False, True]
    ran = [True, False]

    for T, al, ad, r in product(Tarr, alg, add_avoid, ran):
        if ran:
            for i in np.random.choice(fit_idx, size=rep, replace=True):
                yield path, i0, j0, avoid, i, N, T, ad, al
        else:
            for i in range(rep):
                yield path, i0, j0, avoid, -1, N, T, ad, al


def run_evolutionary_params(path, i0=4, j0=6, rep=1000, N=2000):
    fmin = np.load(path)
    p_i = path.stem.split('_')[1]
    seq_str = np.load(f"/data/Jmcbride/RecProt/Data/Sequences/seq_str_{p_i}.npy")
    avoid = np.arange(1, 12, 2)
    fit_idx = np.where(get_fitness_map(fmin, i0, avoid) > 0)[0]
    if len(fit_idx) == 0:
        return None
    with Pool(60) as pool:
        res = np.array(list(pool.starmap(evolve_new_fun, generate_2_(path, i0, j0, avoid, fit_idx, N, rep))))
#   return res
    return res.reshape(7,3,2,2,1000,4,2000)


def run_all_paths():
    base_path = [Path("/data/Jmcbride/RecProt/Data/OpenOpenSize"),
                 Path("/data/Jmcbride/RecProt/Data/LocalStruct11"),
                 Path("/data/Jmcbride/RecProt/Data/LocalStruct12")]

    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']

    for i, ID in enumerate(id_list):
        for base in base_path:
            ts = time.time()
            print(base, ID)
            path = base.joinpath(f"fmin_{i+1}.npy")
            out_path = base.joinpath(f"func_evolve_res_{i+1}.npy")
            if out_path.exists():
                continue
            res = run_evolutionary_params(path)
            if isinstance(res, type(None)):
                continue
            np.save(out_path, res)
            print(f"Time taken: {(time.time()-ts)/60} minutes\n")
    


def mutate(s):
    l = list(s)
    i = np.random.randint(len(s))
    l[i] = {'1':'2', '2':'1'}[l[i]]
    return ''.join(l)


def to_avoid(n, tot):
    return np.random.choice(range(tot), size=n, replace=False)
    

### Multiple fitness functions were tested;
### most give similar results.
### The one used in the paper is "tsvi"
def get_fitness_map(ftot, target, avoid, alg='loose', cut=2, w=1):
    fmin = np.min([np.zeros(ftot.shape), ftot], axis=0)
    if alg == 'loose':
        return np.exp(-ftot[:, target]) - np.sum(np.exp(-ftot[:, avoid]), axis=1)

    elif alg == 'looser':
        f = np.exp(-fmin[:, target]) - np.sum(np.exp(-fmin[:, avoid]), axis=1)
        return np.max([f, np.zeros(f.shape)], axis=0)

    elif alg == 'strict':
        fit =  np.exp(-fmin[:, target])
        favoid = np.min(fmin[:, avoid], axis=1)
#       idx = (favoid <= -cut) | (fmin[:,target] > -cut)
        idx = (favoid <= 0) | (fmin[:,target] > -cut)
        print(f"{np.sum(idx)} sequences removed")
        fit[idx] = 0
        return fit

    elif alg == 'new':
        fit =  np.zeros(fmin.shape[0], float)
        favoid = np.min(fmin[:, avoid], axis=1)
        idx = (fmin[:, target] - favoid < -cut) & (fmin[:,target] < 0)
        fit[idx] = 1
        return fit

    elif alg == 'suggested':
        return 1 / (1 + np.exp(ftot[:, target]) + np.sum(np.exp(ftot[:,[target]] - ftot[:, avoid]), axis=1))

    elif alg == 'tsvi':
        ec = np.exp(-ftot[:, target])
        enc = w * np.sum(np.exp(-ftot[:, avoid]), axis=1)
        return (ec - enc) / (1 + ec + enc)


### Function used instead of cdist;
### this saves memory since cdist cannot work with integers
def get_seqdist_int(S):
    xi = np.array([S]*S.shape[0])
    yi = np.array([[s]*S.shape[0] for s in S])
    return np.abs(xi-yi).sum(axis=2)


def calculate_traversability(seq, fit):
    seq = seq[fit>0].astype("int16")
    print(f"{len(seq)} sequences are fit")
    if len(seq) == 0:
        return np.nan, np.nan, np.nan
    elif len(seq) > 80000:
        return np.nan, np.nan, np.nan

#   sdist = cdist(seq, seq, metric='cityblock')
    sdist = get_seqdist_int(seq)
    clustering = DBSCAN(eps=1, min_samples=1, metric='precomputed').fit(sdist)
    cl = clustering.labels_
    max_traverse = []
    neighbours = np.sum(sdist==1, axis=0).mean()
    for i in np.unique(cl):
        idx = np.where(cl==i)[0]
        for j in idx:
            max_traverse.append(np.max(sdist[j, idx]))
    return np.mean(max_traverse), neighbours, np.unique(cl).size


def get_nearest_noncognate(j, nn, jmax=13, j0=6):
    if nn == 1:
        return [j - nn if j <= j0 else j + nn]
    elif nn == 2:
        return [j - nn if j <= j0 + 1 else j + nn]
    elif nn == 3:
        return [j - nn if j <= j0 + 1 else j + nn]
    

def evolvability_analyses(seq, ftot, cut=2, nn=1, nearest_only=False, fit_cut=0.95):
    alg = 'tsvi'

    targets = np.arange(0, 13, 1)
    avoid_set = [np.arange(1, 12, 2), np.arange(0, 14, 2)]

    frac_fit = np.zeros((targets.size), float)
    traverse = np.zeros((targets.size), float)
    neighbour = np.zeros((targets.size), float)
    nminima = np.zeros((targets.size), float)

    for j, t in enumerate(targets):
        ts = time.time()
        print(f"{alg}\t{j}: {t}")
        if nearest_only:
            # Run only on a single non-cognate, which is nn away from j
            if (j - nn < 0) | (j + nn >= len(targets)):
                continue
            avoid = get_nearest_noncognate(j, nn)
        else:
            # Run on multiple non-cognates, which are within nn of j
            avoid = [k for k in targets if 0 < abs(k-t) <= nn]
        # Get fitness (0-1) and take away 0.5, so that
        # only sequences with fitness > 0.5 are 'fit'
        fit = get_fitness_map(ftot, t, avoid, alg, cut) - fit_cut
        frac_fit[j] = np.sum(fit>0) / fit.size
        out = calculate_traversability(seq, fit)
        traverse[j] = out[0] / len(seq[0])
        neighbour[j] = out[1] / len(seq[1])
        nminima[j] = out[2]
        print(f"Time taken: {(time.time()-ts)/60} minutes\n")
    return frac_fit, traverse, neighbour, nminima


def evaluate_parameter_memory_burden(ftot, cut=2, nearest_only=False, nn=1):
    fmin = np.min([np.zeros(ftot.shape), ftot], axis=0)
    fit_type = ['loose']
    targets = np.arange(0, 13, 1)
    avoid_set = [np.arange(1, 12, 2), np.arange(0, 14, 2)]
    nn_list = np.arange(1, 7)
    nmax = np.zeros((len(fit_type), nn_list.size), float)
    for i, alg in enumerate(fit_type):
        for j, t in enumerate(targets):
            for l, nn in enumerate(nn_list):
#               print(f"{i}: {alg}\t{j}: {t}")
                if nearest_only:
                    avoid = [k for k in targets if 0 < abs(k-t) <= nn]
                else:
                    avoid = avoid_set[j%2]
                fit = get_fitness_map(fmin, t, avoid, alg, cut)
                N = np.sum(fit>0)
                nmax[i,l] = max(nmax[i,l], N)
#               print(alg, t, nn, N)
    return nmax
    

def strict_fit_vs_cut(seq, ftot):
    fmin = np.min([np.zeros(ftot.shape), ftot], axis=0)
    targets = np.arange(0, 14, 2)
    avoid = np.arange(1, 12, 2)
    cut_arr = np.arange(0.5, 10.5, 0.5)
    frac_fit = np.zeros((len(cut_arr), targets.size), float)
    for i, cut in enumerate(cut_arr):
        for j, t in enumerate(targets):
            fit = get_fitness_map(fmin, t, avoid, 'strict', cut)
            frac_fit[i,j] = np.sum(fit>0) / fit.size
    return frac_fit
        

def gen_inp():
#   base_list = ["LocalStruct11", "LocalStruct12", "OpenOpenSize"]
#   nhard = range(1,4)
#   nn_only = [0, 1]
#   id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
#   jsize = range(1, 8)
#   fit_cut = [-0.5, -0.25, 0, 0.25, 0.5]

    base_list = ["LocalStruct11", "LocalStruct12", "OpenOpenSize"]
    nhard = range(1,4)
    nn_only = [0]
    id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-']
    jsize = range(1, 6)
    fit_cut = [0]

    for base, nn, no, i, fc in product(base_list, nhard, nn_only, jsize, fit_cut):
        path = Path(f'../Data/EvolveSim2/new_{base}_{nn}_{i}_{no}_{fc}.pickle')
        path_seq = Path(f'/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/{id_list[i-1]}/Template/inputs/prot_seq.dat')

        yield [path, path_seq, base, nn, no, i, fc]
        

def run_one(inputs):
    path, path_seq, base, nn, no, i, fc = inputs
#   if path.exists():
#       return
    ts = time.time()
    print(path)
    seq = np.loadtxt(path_seq)
    energy = np.load(f"../Data/{base}/ftot_{i}.npy")
    out = evolvability_analyses(seq, energy, nn=nn, nearest_only=no, fit_cut=fc)
    pickle.dump(out, open(path, 'wb'))
    print(f"{path} time: {(time.time()-ts)/60} minutes\n")


def run_all():
    with Pool(5) as pool:
#       _ = pool.starmap(run_one, gen_inp())
#       _ = pool.starmap(run_one, list(gen_inp())[::-1])
        _ = list(pool.imap_unordered(run_one, gen_inp()))



if __name__ == "__main__":

    run_all()

#   nn = 4 - int(sys.argv[1])
#   nn = int(sys.argv[1])
#   base = sys.argv[2]
#   id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-', 'w3h1', 'w2h3']
#   id_list = ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-']#, 'w3h1', 'w2h3']
#   s_list = [np.loadtxt(f'/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/{ID}/Template/inputs/prot_seq.dat') for ID in id_list]
#   for base in ['LocalStruct11', 'LocalStruct12', 'OpenOpenSize']:
#       print(f"Runnning {base}, {nn}")

#   for i in range(7):
#       for no in [0, 1]:
#   for no in [0, 1]:
#       for i in range(7):
#           for fc in [-0.5, -0.25, 0, 0.25, 0.5]:
#               print(f"Runnning {base}, {nn}, {i}, {no}, {fc}")
#               path = Path(f'../Data/EvolveSim2/new_{base}_{nn}_{i+1}_{no}_{fc}.pickle')
#               if path.exists():
#                   continue
#               seq = np.loadtxt(f'/molahome/jmcbride/RecProt/Results/OpenOpenSize/Iter02/{id_list[i]}/Template/inputs/prot_seq.dat')
#               energy = np.load(f"../Data/{base}/ftot_{i+1}.npy")
#               out = evolvability_analyses(seq, energy, nn=nn, nearest_only=no, fit_cut=fc)
#               pickle.dump(out, open(path, 'wb'))
#       
        
    



