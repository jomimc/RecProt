from itertools import product
import os

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix


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



def find_triangles(adj):
    tri = set()
    for i in range(len(adj)):
        i0 = np.where(adj[i])[0]
        for j in i0:
            j0 = np.where(adj[j])[0]
            for k in j0:
                if k != i and k in i0:
                    tri.add(",".join([str(x) for x in sorted([i,j,k])]))
    return np.array([t.split(',') for t in tri], int)


def wedge_product(tri):
    x, y = tri.T
    return np.sign((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))


def bound_dechem(seq, lig1, lig2, emat):
    e1 = np.sum([emat[int(s)-1,int(l)-1] for s, l in zip(seq, lig1)])
    e2 = np.sum([emat[int(s)-1,int(l)-1] for s, l in zip(seq, lig2)])
    return e1 - e2


def bound_dechem_df(df):
    out = []
    for row in df.itertuples():
        out.append(bound_dechem(row.pseq, row.lseq1, row.lseq2, row.emat))
    return out


def get_plot_bins():
    bins = {'ent': np.linspace(4, 11, 50),
            'edef1': np.linspace(0, 10, 51),
            'edef2': np.linspace(0, 10, 51),
            'detot': np.linspace(-10, 10, 81),
            'dedef': np.linspace(-10, 10, 81),
            'dechem': np.linspace(-10, 10, 81),
            'ftot1': np.linspace(-25, 25, 101),
            'ftot2': np.linspace(-25, 25, 101),
            'meanK': np.arange(20, 500, 20),
            'meanE': np.arange(6, 38, 3)}

    combos = [('meanK', 'meanE'),
              ('ftot1', 'ftot2'),
              ('dedef', 'dechem'),
              ('dedef', 'detot'),
              ('edef1', 'detot'),
              ('edef2', 'detot'),
              ('dechem', 'detot')]

    return bins, combos


def eval_input(inp):
    try:
        return eval(inp)
    except:
        return str(inp)


def get_mindist(df):
    pseq = np.array(list(df.pseq.values))
    dist = cdist(pseq, pseq, metric='cityblock')
    return dist.min(axis=0)


def get_neighbours(df):
    pseq = np.array(list(df.pseq.values))
    dist = cdist(pseq, pseq, metric='cityblock')
    np.fill_diagonal(dist, 100)
    return np.sum(dist==1, axis=0)


def update_energy_seq(ener, ang, seq, emat):
    theta = np.pi/3
    for i, s in enumerate(seq):
        e = np.array([emat[1, int(s[i])-1] for i in range(3)])
        idx = ang < theta
        ener[:,i,idx] = np.min([ener[:,i,idx], np.zeros((ener.shape[0], idx.sum())) - max(e[0], e[2])], axis=0)
        idx = ang == theta
        ener[:,i,idx] = np.min([ener[:,i,idx], np.zeros((ener.shape[0], idx.sum())) - np.sum(e)], axis=0)
        idx = ang > theta
        ener[:,i,idx] = np.min([ener[:,i,idx], np.zeros((ener.shape[0], idx.sum())) - max([e[:2].sum(), e[1:].sum()])], axis=0)
    return ener


def get_rec1(e, di=1, dg=2):
    e1 = e[:,di:]
    e2 = e[:,:-di]
    g = np.abs(e1 - e2)
    return (g > dg) & ( ((e1 < 0)&(e2 > 0)) | ((e2 < 0)&(e1 > 0)))


def get_rec2(e, di=1, dg=2):
    e1 = e[:,di:]
    e2 = e[:,:-di]
    g = np.abs(e1 - e2)
    return (g > dg) & ((e1 < 0) | (e2 < 0))


def evaluate_phenotype(e, di=1, dg=2):
    thetaA = np.arange(30, 90, 5) * np.pi/180
    g = e[:,di:] - e[:,:-di]
    r1 = get_rec1(e, di, dg)
    r2 = get_rec2(e, di, dg)

    meang_r1 = [np.mean(np.abs(g[r1[:,i],i])) if np.sum(r1[:,i]) else 0 for i in range(g.shape[1])]
    meang_r2 = [np.mean(np.abs(g[r2[:,i],i])) if np.sum(r2[:,i]) else 0 for i in range(g.shape[1])]

    maxg_r1 = [np.max(np.abs(g[r1[:,i],i])) if np.sum(r1[:,i]) else 0 for i in range(g.shape[1])]
    maxg_r2 = [np.max(np.abs(g[r2[:,i],i])) if np.sum(r2[:,i]) else 0 for i in range(g.shape[1])]

    tot_r1 = r1.sum(axis=0)
    tot_r2 = r2.sum(axis=0)

    frac_r1 = r1.mean(axis=0)
    frac_r2 = r2.mean(axis=0)

    return meang_r1, meang_r2, maxg_r1, maxg_r2, tot_r1, tot_r2, frac_r1, frac_r2




