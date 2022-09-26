from itertools import permutations
import os
from pathlib import Path
import pickle

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from palettable.scientific.diverging import Cork_3, Cork_20, Vik_3, Vik_20
from palettable.scientific.sequential import Devon_20, Bilbao_20, GrayC_20, Acton_16
from palettable.colorbrewer.sequential import Blues_9, YlOrBr_9, Greens_9
import pandas as pd
import seaborn as sns
import trimesh
from vapory import *

import create_inputs as CI
import gene2struct as GS
import rec_io
import rec_utils as utils


def create_protein_object(xy, adj, col='', z=0, radius=0.2, ec=0.2, seq=''):
    edge_col = tuple([ec]*3)
    rad_key = {2: radius*0.05, 3: radius*0.1, 4: radius*0.2}

    objects =[]

    for i, (x, y) in enumerate(xy):
        if isinstance(col, str):
            objects.append(Sphere([x, y, z], radius, Texture(Pigment('color',edge_col)), "no_shadow"))
        else:
            objects.append(Sphere([x, y, z], radius, Texture(Pigment('color',col[i])), "no_shadow"))

    for i, j in zip(*np.where(adj)):
        e1 = [xy[i,0], xy[i,1], 0]
        e2 = [xy[j,0], xy[j,1], 0]
        if len(seq):
            key = int(seq[i]) + int(seq[j])
            r = rad_key[key]
            objects.append(Cylinder(e1, e2, r, Texture(Pigment('color',edge_col)), "no_shadow"))
        else:
            objects.append(Cylinder(e1, e2, radius*.2, Texture(Pigment('color',edge_col)), "no_shadow"))

    return objects


def create_ligand_object(xy, z=0, radius=0.17):
    lig_col = (0.8, 0.8, 0.8)
    edge_col = (0.4, 0.4, 0.4)

    objects =[]

    for i, (x, y) in enumerate(xy):
        objects.append(Sphere([x, y, z], radius, Texture(Pigment('color',lig_col)), "no_shadow"))

    for i in [0,1]:
        j = i + 1
        e1 = [xy[i,0], xy[i,1], 0]
        e2 = [xy[j,0], xy[j,1], 0]
        objects.append(Cylinder(e1, e2, radius*.2, Texture(Pigment('color',edge_col)), "no_shadow"))

    return objects

def light_combos(lvl=0.8):
    lights = [LightSource([ 100,  100,  100], 'color', [lvl, lvl, lvl]),
              LightSource([-100,  100,  100], 'color', [lvl, lvl, lvl]),
              LightSource([ 100, -100,  100], 'color', [lvl, lvl, lvl]),
              LightSource([ 100,  100, -100], 'color', [lvl, lvl, lvl]),
              LightSource([-100,  100, -100], 'color', [lvl, lvl, lvl]),
              LightSource([-100, -100,  100], 'color', [lvl, lvl, lvl]),
              LightSource([ 100, -100, -100], 'color', [lvl, lvl, lvl]),
              LightSource([-100, -100, -100], 'color', [lvl, lvl, lvl])]
    for n in range(1, len(lights)+1):
        for mut in permutations(lights, n):
            yield mut


def create_scene(path, objects, target, source, lights=''):

    camera = Camera('perspective', 'location', source, 'look_at', target)

    if isinstance(lights, str):
        lvl = 0.8
#       lights = [LightSource([ 100,  100,  100], 'color', [lvl, lvl, lvl]),
#                 LightSource([-100,  100,  100], 'color', [lvl, lvl, lvl]),
#                 LightSource([ 100, -100,  100], 'color', [lvl, lvl, lvl]),
#                 LightSource([ 100,  100, -100], 'color', [lvl, lvl, lvl]),
#                 LightSource([-100,  100, -100], 'color', [lvl, lvl, lvl]),
#                 LightSource([-100, -100,  100], 'color', [lvl, lvl, lvl]),
#                 LightSource([ 100, -100, -100], 'color', [lvl, lvl, lvl]),
#                 LightSource([-100, -100, -100], 'color', [lvl, lvl, lvl])]
        lights = [LightSource([ 0,  -20,  0], 'color', [lvl, lvl, lvl]),
                  LightSource([0, 0, -10], 'color', [lvl, lvl, lvl]),
                  LightSource([-10, -10, -10], 'color', [lvl, lvl, lvl])]

    scene = Scene(camera=camera, objects=lights + objects)
    scene.render(path, width=1600, height=1600, antialiasing=0.01, quality=10, output_alpha=True)


def render_protein_feat(path, p_xy, adj, seq, ligands, camera=[0, 0, -6], cmap=Vik_20.mpl_colormap, target=''):
    if isinstance(seq, str):
        col = ''
    elif isinstance(seq, list):
        col = seq
    else:
        col = [cmap(s-1) for s in seq]

    obj = create_protein_object(p_xy, adj, col)
    for l_xy in ligands:
        obj.extend(create_ligand_object(l_xy))
    if isinstance(target, str):
        target = [p_xy[:,0].mean(), p_xy[:,1].mean(), 0]
    create_scene(path, obj, target, camera)


def save_colorbar():
    n_col = 20
    fig, ax = plt.subplots(figsize=(2*n_col,2))
    PC.my_palplot(Vik_20.mpl_colors[::-1], ax=ax)
    ax.set_xticks([-0.5, n_col-0.5])
    x.set_xticklabels([" "*24 + "low / negative", "high / positive"+" "*24], fontsize=72)
    fig.savefig(f'../Figures/prot_colorbar.png', bbox_inches='tight')
    fig.savefig(f'../Figures/prot_colorbar.svg', bbox_inches='tight')



def protein_shapes():
    shapes = CI.protein_shapes()
    for shape in shapes:
        path = Path(f"../../Figures/Drawings/shape_{shape['name']}.png")
        render_protein_feat(path, shape['xy'], shape['adj'], "", [])
        


def create_images():
    c0 = np.loadtxt('../testing/inputs/prot_xy.dat')
    cl = np.loadtxt('../testing/inputs/ligands.dat')[:,4:]
    seq = np.loadtxt('../testing/inputs/prot_seq.dat')
    adj = utils.break_links(utils.get_adj_from_coord(c0))


    seq_idx = np.random.choice(range(len(seq)), size=4)

    for i, c in enumerate(cl[:-3]):
        c1 = c.reshape(2,3).T
        c1[:,1] = c1[:,1] - 0.5
        c1[:,0] = c1[:,0] - 1.3
        c2 = cl[i+3].copy().reshape(2,3).T
        c2[:,1] = c2[:,1] - 0.5
        c2[:,0] = c2[:,0] + 1.3
        print(i, c)
        path = Path(f'../../Figures/Drawings/noseq_widepair_{i}.png')
        render_protein_feat(path, c0, adj, "", [c1, c2])
#       for j, k in enumerate(seq_idx):
#           print(j, k)
#           path = Path(f'../../Figures/Drawings/basic_{i}_{j}.png')
#           render_protein_feat(path, c0, adj, seq[k], [c])



def fig1():
    c0 = np.loadtxt('../testing/inputs/prot_xy.dat')
    adj = utils.break_links(utils.get_adj_from_coord(c0))

    angles = np.array([80, 60, 40]) * np.pi / 180
    dx = [-2, 0, 2]
    dy = [-0.2 - np.sin(x/180*np.pi) for x in [20, 30, 40]]
    lig = []
    for i, a in enumerate(angles):
        c = CI.basic_ligand(a)
        c[:,0] = c[:,0] + dx[i]
        c[:,1] = c[:,1] + dy[i]
        lig.append(c)

    camera = [0, -0.5, -6.5]
    target = [0, -0.25, 0]
    path = Path(f'/home/johnmcbride/LaTEX/RecProtPaper/Figures/fig1.png')

    b_idx = [1, 6, 2]
    cmap = Acton_16.mpl_colormap
    col_seq = [cmap(0.5)[:3] if i in b_idx else (0.4, 0.4, 0.4) for i in range(len(adj))]
    render_protein_feat(path, c0, adj, col_seq, lig, camera=camera, target=target)
        

def fig3():
    camera = [0, -0.5, -2.5]
    target = [0, -0.25, 0]
    angles = np.array([45, 40]) * np.pi / 180
    for i, a in enumerate(angles):
        c = CI.basic_ligand(a)
        path = Path(f'../../Figures/Drawings/fig3_{i}.png')
        obj = create_ligand_object(c)
        create_scene(path, obj, target, camera)


def create_custom_protein(xy, adj, seq, col, z=0, radius=0.2, ec=0.6):
    edge_col = tuple([ec]*3)

    objects =[]

    for i, (x, y) in enumerate(xy):
        s1 = Sphere([x, y, z], radius, Texture(Pigment('color',col[i])), "no_shadow")
        b1 = Box([x, y + radius, z + radius], [x + radius, y - radius, z - radius])
        objects.append(Difference(s1, b1))

        if i in [1, 2, 6]:
            i2 = {1:-3, 2:-1, 6:-2}[i]
            s2 = Sphere([x, y, z], radius, Texture(Pigment('color',col[i2])), "no_shadow")
        else:
            s2 = Sphere([x, y, z], radius, Texture(Pigment('color', edge_col)), "no_shadow")
        b2 = Box([x, y + radius, z + radius], [x - radius, y - radius, z - radius])
        objects.append(Difference(s2, b2))

    col_key = {2: (0.8, 0.8, 0.8), 3: (0.6, 0.6, 0.6), 4: (0.4, 0.4, 0.4)}
    rad_key = {2: radius*0.05, 3: radius*0.1, 4: radius*0.2}

    for i, j in zip(*np.where(adj)):
        e1 = [xy[i,0], xy[i,1], 0]
        e2 = [xy[j,0], xy[j,1], 0]

        key = int(seq[i]) + int(seq[j])
        c = col_key[key]
        r = rad_key[key]
#       objects.append(Cylinder(e1, e2, r, Texture(Pigment('color', c)), "no_shadow"))
        objects.append(Cylinder(e1, e2, r, Texture(Pigment('color', edge_col)), "no_shadow"))

    return objects


def fig3a():
    camera = [0, -0.5, -6.5]
    target = [0, -0.25, 0]
    c0 = np.loadtxt('../testing/inputs/prot_xy.dat')
    adj = utils.break_links(utils.get_adj_from_coord(c0))

    lights = list(light_combos())
#   idx = np.random.choice(range(len(lights)), replace=False, size=50)
#   for i in idx:
#       path = Path(f"../../Figures/Drawings/tmp/{i:06d}.png")
#       cx = {'1':0.2, '2':0.8}
#       seq = "2112222122212222"
#       col = [list(Vik_20.mpl_colormap(cx[s])) for s in seq]
#       col = [c[:3] + [0.50] for c in col]
#       obj = create_custom_protein(c0, adj, seq, col)
#       create_scene(path, obj, target, camera, lights=list(lights[i]))

    path = Path(f"../../Figures/Drawings/fig3_prot.png")
#   cx = {'1':list(Vik_20.mpl_colormap(0.2)), '2':list(Vik_20.mpl_colormap(0.8)),
    cx = {'1':list(sns.color_palette('dark')[0]) +[1.0],
          '2':list(sns.color_palette('dark')[1]) +[1.0],
          '3':list(sns.color_palette('deep')[2]) +[1.0]}
    seq = "2112222122212333"
    col = [cx[s] for s in seq]
    col = [c[:3] + [0.50] for c in col]
    obj = create_custom_protein(c0, adj, seq, col, ec=0.7)
    l = list(lights[31761])
    create_scene(path, obj, target, camera, lights=l)


def fig3b():
    camera = [0, -0.5, -6.5]
    target = [0, -0.25, 0]
    c0 = np.loadtxt('../testing/inputs/prot_xy.dat')
    adj = utils.break_links(utils.get_adj_from_coord(c0))

    lights = list(light_combos())

    path = Path(f"../../Figures/Drawings/fig3_prot2.png")
#   cx = {'1':list(Vik_20.mpl_colormap(0.2)), '2':list(Vik_20.mpl_colormap(0.8)),
    cx = {'1':list(sns.color_palette('dark')[0]) +[1.0],
          '2':list(sns.color_palette('dark')[1]) +[1.0],
          '3':list(sns.color_palette('deep')[2]) +[1.0]}
    seq = "2112222122212"
    col = [cx[s][:3] + [0.00] for s in seq]
#   obj = create_custom_protein(c0, adj, seq, col, ec=0.7)
    obj = create_protein_object(c0, adj, col, ec=0.7, seq=seq)
    l = list(lights[31761])
    create_scene(path, obj, target, camera, lights=l)


def fig4():
    camera = [0, -0.5, -8.5]
    target = [0, -0.25, 0]
    shapes = CI.protein_shapes()
    j = 0
    col = np.array(sns.color_palette())[[0,1,3,2,4]]
    for i, shape in enumerate(shapes):
        if shape['name'] in ['w1h1', 'w2h1', 'w2h2-', 'w2h2', 'w2h3-']:
            path = Path(f"../../Figures/Drawings/fig4_{shape['name']}.png")
            c = [col[j]]*len(shape['xy'])
            obj = create_protein_object(shape['xy'], shape['adj'], c)
            create_scene(path, obj, target, camera)
            j += 1


def fig5():
    camera = [0, 0.0, -1.0]
    target = [0, 0.0, 0]
    c0 = np.loadtxt('../testing/inputs/prot_xy.dat')
    adj = utils.break_links(utils.get_adj_from_coord(c0))

    cx = [list(sns.color_palette('dark')[0]),
          list(sns.color_palette('dark')[1]),
          list(sns.color_palette('dark')[6]),
          list(sns.color_palette('deep')[2])]
    lights = list(list(light_combos())[31761])
#   lights = list(list(light_combos())[-1])
    for i in [0,1,2,3]:
        path = Path(f"../../Figures/Drawings/fig5_{i}.png")
        seq = str(i+1)
        col = [cx[i]]
#       obj = create_custom_protein([[0,0]], [], [seq], col, ec=0.7)
        obj = create_protein_object([[0,0]], [], col, ec=0.7)
        create_scene(path, obj, target, camera, lights=lights)

    path = Path(f"../../Figures/Drawings/fig5_prot.pdf")
    cx = {1:0.1, 2:0.9}
    seq = np.array(list("2112222122212")).astype(int)

    x = 0.3
    tab = np.array([[-x, 0], [0, x]])
    c1 = GS.get_new_coord(seq[:13]-1, adj, c0, tab)

    fig, ax = plt.subplots(figsize=(8,8))
    for i in range(13):
#       ax.plot([c0[i,0]], [c0[i,1]], 'o', c=Vik_20.mpl_colormap(cx[seq[i]]), fillstyle='left', mec=(0.5, 0.5, 0.5), mew=1.5, ms=15)
#       ax.plot([c0[i,0]], [c0[i,1]], 'o', c=Vik_20.mpl_colormap(cx[seq[i]]), fillstyle='none', mec=(0.5, 0.5, 0.5), mew=1.5, ms=10)
        ax.plot([c0[i,0]], [c0[i,1]], 'o', c='k', fillstyle='none', ms=16)
#       ax.plot([c1[i,0]], [c1[i,1]], 'o', c=Vik_20.mpl_colormap(cx[seq[i]]), fillstyle='left', mec=(0.5, 0.5, 0.5), mew=1.5, ms=20)
#       ax.plot([c0[i,0], c1[i,0]], [c0[i,1], c1[i,1]], '-k', lw=2)
#       ax.quiver(c0[i,0], c0[i,1], c1[i,0] - c0[i,0], c1[i,1] - c0[i,1])

    for d in ['left', 'right', 'top', 'bottom']:
        ax.spines[d].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-1.5, 3.5)

    fig.savefig(path, bbox_inches='tight', transparent=True)



if __name__ == "__main__":

#   fig3b()
    fig5()



