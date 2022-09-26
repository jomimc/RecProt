from itertools import product
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

import rec_utils as utils


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


def getK(A, kmat, seq):
    seq = np.array(list(seq), int) - 1
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if i > j])
    Nb = len(bonds)
    K = np.zeros((Nb, Nb), float)
    np.fill_diagonal(K, [kmat[seq[i],seq[j]] for i, j in zip(*bonds.T)])
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


### Add bond at binding site to stiffen binding site;
### Change the strength of the bonds at the binding site to "kstiff"
def add_bonds(K, D, A, pc, idx0=[1,2,6], kstiff=10):
    nij = getnij(pc)

    Nbold = len(K)
    Nb = len(K) + 1
    Knew = np.zeros((Nb, Nb), float)
    Knew[:Nbold,:Nbold] = K.copy()

    idx = idx0 + [Nb-1]
    for i in idx:
        Knew[i,i] = kstiff

    A = A.copy()
    A[idx[0],idx[2]] = 1
    A[idx[2],idx[0]] = 1

    grad = get_grad(A)

    return A, grad, Knew, getD(grad, nij)


def get_H(K, D):
    return D.T.dot(K).dot(D)


### MAIN Function
### for calculating the entropy change upon binding
def entropy_diff(ecoord, adj, kmat, seq, idx0=[1,6,2], kstiff=10):
    K = getK(adj, kmat, seq)
    nij = getnij(ecoord)
    grad = get_grad(adj)
    D = getD(grad, nij)
    H = get_H(K, D)

    val, vec = np.linalg.eig(H)
    # Ignore the trivial zero-energy modes (translation, rotation)
    p1 = np.product([v.real for v in val if v > 10**-5])

    An, Gn, Kn, Dn = add_bonds(K, D, adj, ecoord, idx0=idx0, kstiff=kstiff)
    Hn = get_H(Kn, Dn)
    valn, vecn = np.linalg.eig(Hn)
    # Ignore the trivial zero-energy modes (translation, rotation)
    p2 = np.product([v.real for v in valn if v > 10**-5])
    return np.log(p2/p1)/2


def get_kmat(top):
    lo = top * 0.1
    mid = top * 0.5
    return np.array([[lo, mid], [mid, top]])


def fromKtoeigval(adj, ecoord, seq, kstiff=10):
    smax = 10.**np.linspace(0, 2, 30)
    meaneig = []
    ent = []
    ratio = []
    for s in smax:
        kmat = get_kmat(s)
        K = getK(adj, kmat, seq)
        nij = getnij(ecoord)
        grad = get_grad(adj)
        D = getD(grad, nij)
        H = get_H(K, D)

        val, vec = np.linalg.eig(H)
        meaneig.append(np.mean(val[:-3].real))
        p1 = np.product([v.real for v in val if v > 10**-6])
        ent.append(-0.5 * np.log(p1) + (1 + np.log(np.pi * 2)) * (len(val) - 3) / 2)

        An, Gn, Kn, Dn = add_bonds(K, D, adj, ecoord, stiff=kstiff)
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
        adj = utils.break_links(utils.get_adj_from_coord(en))
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



