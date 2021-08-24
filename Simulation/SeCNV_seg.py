import sys
import numpy as np
import time
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
from scipy.spatial.distance import pdist, squareform
from numpy import exp
import KR_norm_juicer

file_name = sys.argv[1]

def load_matrix(file_name):
    return np.loadtxt(file_name)

class DP_process:
    def __init__(self,input_matrix,K_max):
        self.matrix = input_matrix
        self.n = len(input_matrix)
        self.max_k = K_max
        self.dp_table = np.zeros((self.n, self.max_k))
        self.index_table = np.zeros((self.n, self.max_k))
        self.pre_table = np.zeros((self.n, self.n))
        for start in range(0, self.n):
            for end in range(start, self.n):
                current_matrix = self.matrix[start:end + 1, start:end + 1]
                self.pre_table[start][end] = 2 * np.sum(np.triu(current_matrix,1))
                self.pre_table[end][start] = np.sum(self.matrix[0:start, start:end + 1]) \
                                             + np.sum(self.matrix[end + 1:self.n, start:end + 1])
        self.sum = self.pre_table[0][self.n - 1]
        self.dlogd_sum = np.zeros(self.n)
        self.dlogd_sum[0] = self.pre_table[0][0] * np.log2(self.pre_table[0][0])
        for bin in range(1, self.n):
            self.dlogd_sum[bin] = self.dlogd_sum[bin - 1] + self.pre_table[bin][bin] * np.log2(self.pre_table[bin][bin])

    def calculate_se(self, g, V_p, V):
        if (V == 0):
            return 0
        else:
            return g / self.sum * np.log2(V_p / V)

    #record the initial state of dynamic programming (K = 0)
    def dp_init(self):
        for i in range(0, self.n):
            if i == 0:
                current_volume = self.pre_table[0][0]
            else:
                current_volume = self.pre_table[0][i] + self.pre_table[i][0]
            intra_se = (current_volume * np.log2(current_volume) - self.dlogd_sum[i])/self.sum
            inter_se = self.calculate_se(self.pre_table[i][0], self.sum, current_volume)
            self.dp_table[i][0] = inter_se + intra_se

    #dynamic process
    def dp_process(self):
        for temp_k in range(1, self.max_k):
            for temp_n in range(0, self.n):
                min_se = np.inf
                min_index = 0
                for i in range(0,temp_n):
                    se_1 =  self.dp_table[i][temp_k - 1]  # T(i,k-1)
                    if i + 1 == temp_n:
                        current_volume = self.pre_table[temp_n][temp_n]
                    else:
                        current_volume = self.pre_table[i+1][temp_n] + self.pre_table[temp_n][i+1]
                    intra_se = (current_volume * np.log2(current_volume) - (self.dlogd_sum[temp_n]-self.dlogd_sum[i])) / self.sum
                    inter_se = self.calculate_se(self.pre_table[temp_n][i + 1], self.sum, current_volume)
                    se_2 = intra_se + inter_se  # H(i+1;n)
                    temp_se = se_1 + se_2
                    if temp_se < min_se:
                        min_se = temp_se
                        min_index = i
                self.dp_table[temp_n][temp_k] = min_se
                self.index_table[temp_n][temp_k] = min_index

    #find the best K
    def find_k(self):
        k_list = self.dp_table[-1].tolist()
        self.k = k_list.index(np.min(k_list)) + 1
        self.min_entropy  = np.min(k_list)


    #find the boundaries:
    def find_boundaries(self):
        self.boundaries = []
        self.boundaries.append(int(self.n - 1))
        for i in range(1, self.k):
            self.boundaries.append(int(self.index_table[self.boundaries[i-1]][self.k-i]))
        self.boundaries.reverse()
        return self.boundaries

    def show_result(self):
        print("The optimal k: {}".format(self.k))
        plt.clf()
        sns.set()
        cmap = sns.cubehelix_palette(gamma=0.8, as_cmap=True)
        sns.heatmap(self.matrix, fmt="d", cmap=cmap)
        #plt.axis('off')
        currentAxis = plt.gca()
        for i in range(self.k):
            if i == 0:
                print("cluster {}:from {} to {}".format(i + 1, 1, self.boundaries[i] + 1))
                rect = patches.Rectangle((0, 0), self.boundaries[i] + 1, self.boundaries[i] + 1,
                                         linewidth=2, edgecolor='royalblue', facecolor='none')
                currentAxis.add_patch(rect)
            else:
                print("cluster {}:from {} to {}".format(i + 1, self.boundaries[i - 1] + 2, self.boundaries[i] + 1))
                rect = patches.Rectangle((self.boundaries[i - 1] + 1, self.boundaries[i - 1] + 1),
                                         self.boundaries[i] - self.boundaries[i - 1], self.boundaries[i] - self.boundaries[i - 1],
                                         linewidth=2, edgecolor='royalblue', facecolor='none')
                currentAxis.add_patch(rect)
        print("Entropy is : {}".format(self.min_entropy))


def make_adj_matrix(matrix, K, sigma):
    if K == "auto_set":
        K = int(0.05 * matrix.shape[1])
    else:
        K = int(K)
    if sigma == "auto_set":
        sigma = np.median(matrix)/2
    else:
        sigma = float(sigma)

    adj_matrix = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i in range(matrix.shape[0]):
        for j in range(i+1, matrix.shape[0]):
            temp = (matrix[i]-matrix[j]) ** 2
            temp.sort()
            temp = temp[-K:]
            adj_matrix[i][j] = np.exp(-np.sum(temp)/sigma**2)
            adj_matrix[j][i] = np.exp(-np.sum(temp)/sigma**2)
    adj_matrix = KR_norm_juicer.KR_norm(adj_matrix)

    return adj_matrix


def main():
    matrix = load_matrix(file_name)
    adj_matrix = make_adj_matrix(matrix, "auto_set", "auto_set")
    dp_process = DP_process(adj_matrix, 20)
    dp_process.dp_init()
    dp_process.dp_process()
    dp_process.find_k()
    boundaries = dp_process.find_boundaries()
    dp_process.show_result()


if __name__ == "__main__":
    main()

