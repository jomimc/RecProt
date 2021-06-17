from itertools import product
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

PATH_BASE = '/home/johnmcbride/Downloads/specificity_john/test'


def add_atoms(ecoord, Nnew):
    N = len(ecoord)
    enew = []
    dx = 1
    dy = 0.5 * 3**0.5
    for i in range(N, N + Nnew):
        j = i % 9
        if j >= 4:
            yi = 2 * (i // 9)
            xi = (j - 4) - 2
        else:
            yi = 2 * (i // 9) - 1
            xi = (j) - 1.5

        enew.append([xi * dx, yi * dy])

    return np.append(ecoord, enew, axis=0)


def get_adj_from_coord(c, d0=1.01):
    dist = distance_matrix(c, c)
    np.fill_diagonal(dist, dist.max())
    return (dist<=d0).astype(int)


def break_links(adj, idx=[[1,2]]):
    for i, j in idx:
        adj[i,j] = 0
        adj[j,i] = 0
    return  adj


def get_equilcoord(base=''):
    x = np.loadtxt(os.path.join(base, 'inputwedge/x13.dat'))
    y = np.loadtxt(os.path.join(base, 'inputwedge/y13.dat'))
    return np.array([x, y]).T


def parse_data(data, ns=''):
    if isinstance(ns, str):
        ns = (2**16, 6, 6)
    iprot, ilig, iang = [data[:,i].astype(int).reshape(ns) for i in range(3)]
    iprot = iprot[:,0,0]
    ilig = ilig[0,:,0]
    iang = iang[0,0,:]
    clig = data[:,3:9].reshape((2**16, 6, 6, 6))
    seq = data[:,9].astype(int).astype(str).reshape(ns)[:,0,0]
    etot, edef, echem = [data[:,i].reshape(ns) for i in range(10, 13)]
    cprot = data[:,13:].reshape(2**16, 6, 6, 26)
    return seq, iprot, ilig, iang, clig, cprot, etot, edef, echem


def prot_df(seq, fmin, edef, echem, ent, gap, lig):
    emin = fmin - ent
    return pd.DataFrame(data={'seq':seq, 'fmin':fmin, 'emin':emin, 'edef':edef,
                              'echem':echem, 'ent':ent, 'gap':gap, 'ligand':lig})


def get_gap(f):
    f = f.ravel()
    f = sorted(f)
    return f[1] - f[0]


def evaluate_data(seq, ftot, etot, edef, echem, ent, adj, kmat, emat, idx):
    Nl = len(idx)*6
    lig_key = {}
    count = 0
    for i in range(6):
        for j in range(6):
            if j in idx:
                lig_key[count] = i*6 + j
                count += 1

    idx = sorted(list(lig_key.values()))

    N = len(ftot)
    ftot = ftot.copy().reshape(N, 36)
    etot = etot.copy().reshape(N, 36)
    edef = edef.copy().reshape(N, 36)
    echem = echem.copy().reshape(N, 36)


    df = pd.DataFrame(data={'seq':seq})
    ekey = {'1':emat[0], '2':emat[1]}
    df['meanE'] = df.seq.apply(lambda x: sum([ekey[y] for y in x[-3:]]))
    df['meanK'] = df.seq.apply(lambda x: np.sum(getK(adj, kmat, x[:-3])))

    df['gap'] = np.array([get_gap(f[idx]) for f in ftot])
    df['ent'] = ent

    best_lig = np.array([[lig_key[j] for j in np.argsort(ftot[i,idx])] for i in range(len(ftot))])

    fmin = np.array([e[i] for e, i in zip(ftot, best_lig)])
    emin = np.array([e[i] for e, i in zip(etot, best_lig)])
    chemmin = np.array([e[i] for e, i in zip(echem, best_lig)])
    defmin = np.array([e[i] for e, i in zip(edef, best_lig)])

    lbls = ['lig', 'fmin', 'emin', 'echem', 'edef']
    for i in range(2):
        for l, d in zip(lbls, [best_lig, fmin, emin, chemmin, defmin]):
            df[f"{l}_{i+1}"] = d[:,i]

    return df



def results_summary(path):
    data = np.loadtxt(os.path.join(PATH_BASE, path, "output.dat"), skiprows=2)
    ks, kw, kww, els, elw = [float(x) for x in path.split('_')]
    kmat = np.array([[kww, kw], [kw, ks]])
    emat = [elw, els]
    kstiff = ks * 10.

    seq, iprot, ilig, iang, clig, cprot, etot, edef, echem = parse_data(data)
    ecoord = get_equilcoord(os.path.join(PATH_BASE, path))
    adj = break_links(get_adj_from_coord(ecoord))
    ent = np.array([entropy_diff(ecoord, clig, adj, kmat, s, kstiff=kstiff) for s in seq])
    ftot = etot + ent[:,None, None]

    df1 = evaluate_data(seq, ftot, etot, edef, echem, ent, adj, kmat, emat, [0,1,2,3,4,5])
    df2 = evaluate_data(seq, ftot, etot, edef, echem, ent, adj, kmat, emat, [0,1,2,3,5])

    df1.to_pickle(os.path.join(PATH_BASE, path, "all.pickle"))
    df2.to_pickle(os.path.join(PATH_BASE, path, "mismatch.pickle"))
    out = [df1.lig_1.value_counts(), df2.lig_1.value_counts()]

    pickle.dump(out, open(os.path.join(PATH_BASE, path, "bl_count.pickle", 'wb')))



def get_A(N=13):
    path = f"inputwedge/adj_list{N}.dat"
    Z = np.zeros((N,N), int)
    A = np.zeros((N,N), int)
    for l in open(path, 'r'):
        i, j = [int(x) - 1 for x in l.strip('\n').split()[:2]]
        Z[i,j] = -1
        Z[j,i] = -1
        Z[i,i] += 1
        Z[j,j] += 1
        A[i,j] = 1
        A[j,i] = 1
    return Z, A


def get_grad(A):
    Nb = int(np.sum(A==1)/2)
    grad = np.zeros((Nb, len(A)), int)
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if j > i])
    for k, (i, j) in enumerate(zip(*bonds.T)):
        grad[k,i] = 1
        grad[k,j] = -1
    return grad


def getK(A, smat, seq):
    seq = np.array(list(seq), int) - 1
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if i > j])
    Nb = len(bonds)
    K = np.zeros((Nb, Nb), float)
    np.fill_diagonal(K, [smat[seq[i],seq[j]] for i, j in zip(*bonds.T)])
    return K


def getnij(coord):
    N = len(coord)
    nij = np.zeros((N, N, 2), float)
    for i, j in zip(*np.tril_indices(N, k=-1)):
        rij = coord[i] - coord[j]
        rij = rij / np.linalg.norm(rij)
        nij[i,j] = rij
        nij[j,i] = -rij
    return nij


def getD(grad, nij, d=2):
    Nb = np.sum(grad==1)
    Na = nij.shape[1]
    D = np.zeros((Nb, int(Na*d)), float)
    for i in range(Nb):
        d0 = np.zeros((Na,d), float)
        j, k = np.where(grad[i] != 0)[0]
        d0[j] = grad[i,j] * nij[j,k]
        d0[k] = grad[i,k] * nij[j,k]
        D[i] = d0.reshape(Na*d)
    return D


def add_bonds(K, D, A, pc, lc, kstiff=10):
    nij = getnij(pc)

    Nbold = len(K)
    Nb = len(K) + 1
    Knew = np.zeros((Nb, Nb), float)
    Knew[:Nbold,:Nbold] = K.copy()

    idx = [1,2,6, Nb-1]
    for i in idx:
        Knew[i,i] = kstiff

    A = A.copy()
    A[1,2] = 1
    A[2,1] = 1

    grad = get_grad(A)

    return A, grad, Knew, getD(grad, nij)


def get_H(K, D):
    return D.T.dot(K).dot(D)


def entropy_diff(ecoord, clig, adj, smat, seq, kstiff=10):

    K = getK(adj, smat, seq)
    nij = getnij(ecoord)
    grad = get_grad(adj)
    D = getD(grad, nij)
    H = get_H(K, D)

    val, vec = np.linalg.eig(H)
#   p1 = np.product(val[:-3].real)
    p1 = np.product([v.real for v in val if v > 10**-5])

    An, Gn, Kn, Dn = add_bonds(K, D, adj, ecoord, clig[0,0,3].reshape(3,2), kstiff=kstiff)
    Hn = get_H(Kn, Dn)
    valn, vecn = np.linalg.eig(Hn)
#   p2 = np.product(valn[:-3].real)
    p2 = np.product([v.real for v in valn if v > 10**-5])
    return np.log(p2/p1)/2


def get_smat(top):
    lo = top * 0.1
    mid = top * 0.5
    return np.array([[lo, mid], [mid, top]])


def fromKtoeigval(adj, ecoord, clig, seq, kstiff=10):
    smax = 10.**np.linspace(0, 2, 30)
    meaneig = []
    ent = []
    ratio = []
    for s in smax:
        smat = get_smat(s)
        K = getK(adj, smat, seq)
        nij = getnij(ecoord)
        grad = get_grad(adj)
        D = getD(grad, nij)
        H = get_H(K, D)

        val, vec = np.linalg.eig(H)
        meaneig.append(np.mean(val[:-3].real))
        p1 = np.product([v.real for v in val if v > 10**-6])
        ent.append(-0.5 * np.log(p1) + (1 + np.log(np.pi * 2)) * (len(val) - 3) / 2)

        An, Gn, Kn, Dn = add_bonds(K, D, adj, ecoord, clig.reshape(3,2), kstiff=kstiff)
        Hn = get_H(Kn, Dn)
        valn, vecn = np.linalg.eig(Hn)
        p2 = np.product([v.real for v in valn if v > 10**-6])
#       p2 = -0.5 * np.log(p2) + np.pi * 2 / (len(val) - 3)
        ratio.append(-0.5 * np.log2(p1 / p2))

    return smax, meaneig, ent, ratio


def vs_size(ecoord, clig, seq, kstiff=10):
    n_arr = np.cumsum([5, 4]*5)
    smax = 10.**np.linspace(0, 2, 30)

    meaneig, ent, ratio = [], [], []
    for n in n_arr:
        en = add_atoms(ecoord, n)
        adj = break_links(get_adj_from_coord(en))
        seq = ''.join([seq[0]]*(len(en)+3))

        smax, me, e, r= fromKtoeigval(adj, en, clig, seq, kstiff=kstiff)
        meaneig.append(me)
        ent.append(e)
        ratio.append(r)

    return n_arr, smax, np.array(meaneig), np.array(ent), np.array(ratio)


def get_quadrant(row, gap_cut=0.2):
    gap, e, f1, f2 = [row[i] for i in [3, 4, 6, 11]]
    return int(gap > gap_cut) + int(f2 > 0) + int(f1 < 0)


def tune_entropy_max_rec(df):
    de_list = np.linspace(-5, 10, 51)
    out = []
    for de in de_list:
        df2 = df.copy()
        df2['ent'] = df2['ent'] + de
        df2['fmin_1'] = df2['fmin_1'] + de
        df2['fmin_2'] = df2['fmin_2'] + de
        out.append(sum(df2.apply(get_quadrant, axis=1)==3))
    return de_list, np.array(out)


def plot_best_rec_pos(df, ec, cprot, clig, gc=0.2):
    for i, c in df.loc[(df.gap>gc)&(df.fmin_2>0), 'lig_1'].value_counts().items():
        idx = (df.gap>gc) & (df.fmin_2>0) & (df.fmin_1<0) & (df.lig_1==i)
        c = sum(idx)
        j = i % 6
        fig, ax = plt.subplots()
        colors = np.array([list(str(x)) for x in df.loc[idx, 'seq']], float)[:,:-3].mean(axis=0)
#       err = np.array([list(str(x)) for x in df.loc[idx, 'seq']], float)[:,:-3].std(axis=0)
        mean_pos = cprot[idx,5,j].mean(axis=0).reshape(13,2)
        ax.plot(*ec.T, 'ok', alpha=0.5)
        for k in range(len(ec)):
            ax.plot([ec[k,0], mean_pos[k,0]], [ec[k,1], mean_pos[k,1]], '-k', alpha=0.5)
        sc = ax.scatter(*mean_pos.T, c=colors)
#       ax.errorbar(*cprot[idx,5,j].mean(axis=0).reshape(13,2).T, yerr=err, linestyle='none')
        ax.plot(*clig[0,5,j].reshape(3,2).T, 's', ms=10, fillstyle='none')
        fig.colorbar(sc)
        ax.set_title(f"{i}: {c}")



