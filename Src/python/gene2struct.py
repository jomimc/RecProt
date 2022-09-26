
import matplotlib.pyplot as plt
import numpy as np

def get_fold_table(dr):
    sigma = np.random.choice(dr, replace=True, size=3)
    return np.array([[sigma[0], sigma[1]], [sigma[1], sigma[2]]])


def get_new_coord(seq, adj, ec, tab):
    nec = ec.copy()
    for i, s in enumerate(seq):
        idx = np.where(adj[i])[0]
        vec = np.mean([tab[seq[i], seq[j]] * (ec[j] - ec[i]) for j in idx], axis=0)
        print(i, s)
        print(idx)
        print(vec)
        nec[i] = nec[i] + vec
    return nec


def plot_alt_coord(ec, adj):
    dr = [-0.2, -0.1, 0, 0.1, 0.2]
    fig, ax = plt.subplots(2,2); ax = ax.reshape(ax.size)
    seq = np.random.randint(2,size=13)
    for i in range(4):
        tab = get_fold_table(dr)
        nec = get_new_coord(seq, adj, ec, tab)
        ax[i].plot(*ec.T, 'ok')
        xc = nec[:,0] - ec[:,0]
        yc = nec[:,1] - ec[:,1]
#       ax[i].quiver(ec[:,0], ec[:,1], xc, yc)
        ax[i].plot(*nec.T, 'o', alpha=0.5)



