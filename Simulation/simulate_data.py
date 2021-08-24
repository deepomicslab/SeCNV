import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
from matplotlib import colors

norm = float(sys.argv[1])
bp = float(sys.argv[2])
gn = float(sys.argv[3])

def make_syn_data(sc_num, bin_num, loc_list, norm_ratio, bp_ratio):
    matrix = np.ones((sc_num, bin_num)) * 2
    norm_num = int(sc_num * norm_ratio)
    for i in range(matrix.shape[0]):
        if i > norm_num:
            p = np.array([0.05, 0.1, 0.06, 0.3, 0.45, 0.02, 0.02])
            cn = np.random.choice([1, 2, 3, 4, 5, 6, 7], p=p.ravel())
            for j in range(matrix.shape[1]):
                if j in loc_list:
                    exchange_pro = random.uniform(0, 1)
                    if exchange_pro > 1 - bp_ratio:
                        p = np.array([0.05, 0.1, 0.06, 0.3, 0.45, 0.02, 0.02])
                        cn = np.random.choice([1, 2, 3, 4, 5, 6, 7], p=p.ravel())
                matrix[i][j] = cn

    return matrix

def draw_cn(matrix, name):
    plt.figure(figsize=(8, 4.5))
    plt.clf()
    sns.set()
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i][j] > 7:
                matrix[i][j] = 7
    cmap = colors.ListedColormap(
        ['lightsteelblue', 'white', 'papayawhip', 'peachpuff', 'sandybrown', 'tomato', 'darkred'])
    sns.heatmap(matrix, fmt="d", cmap=cmap, cbar_kws={"ticks":np.arange(1,8), "boundaries":np.arange(0.5,8.5)})
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.yticks(rotation=0)
    plt.savefig(name, dpi=600)
    plt.close()

def make_read_matrix(cn_matrix, var):
    rn_matrix = np.zeros(cn_matrix.shape)
    for i in range(cn_matrix.shape[0]):
        for j in range(cn_matrix.shape[1]):
            rn_matrix[i][j] = np.random.poisson(lam=cn_matrix[i][j]*200) + np.random.normal(loc=0, scale=var)
    rn_matrix = np.maximum(0, rn_matrix)
    return rn_matrix

def draw_rn(rn_matrix, name):
    plt.figure(figsize=(8, 4.5))
    sns.set()
    cmap = sns.cubehelix_palette(gamma=0.8, as_cmap=True)
    sns.heatmap(rn_matrix, fmt="d", cmap=cmap)
    plt.savefig(name, dpi=600)
    plt.close()

def save_data(rn_matrix, file_name):

    f = open(file_name, "w")
    for i in rn_matrix:
        for j in i:
            f.write(str(j))
            f.write('\t')
        f.write('\n')
    f.close()

cn_matrix = make_syn_data(100, 200, [20, 40, 50, 70, 80, 120, 140, 150, 170, 180], norm, bp)
rn_matrix = make_read_matrix(cn_matrix, gn)
save_data(rn_matrix.T, "data_%s_%s_%s.txt"%(norm, bp, gn))
draw_cn(cn_matrix, "data_%s_%s_%s.png"%(norm, bp, gn))
