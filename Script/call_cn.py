import os
import sys
import numpy as np
import pandas as pd
from numpy import exp
from scipy import stats
import KR_norm_juicer
from sklearn.neighbors import KernelDensity
from scipy import signal, stats
import warnings
warnings.filterwarnings('ignore')

work_path = sys.argv[1]
min_ploidy = float(sys.argv[2])
max_ploidy = float(sys.argv[3])
K = int(sys.argv[4])

def read_matrix(file_name):
    f = open(file_name)
    matrix = []
    chr_name = []
    bin_list = []
    count = 0
    for line in f:
        line = line.strip("\n").rstrip("\t").split("\t")
        count += 1
        if not count == 1:
            matrix.append(line[3:])
            chr_name.append(line[0])
            bin_list.append(line[0]+":"+line[1]+"-"+line[2])
        else:
            sample_list = line
    matrix = np.array(matrix)
    bin_list = np.array(bin_list)
    sample_list = np.array(sample_list[3:])
    matrix = matrix.astype(np.int)
    for i in range(len(sample_list)):
        sample_list[i] = str(sample_list[i]).split("/")[-1]
    chr_name = np.array(chr_name)
    f.close()
    return matrix, chr_name, bin_list, sample_list

def sose(scnp):
    rcnp = np.round(scnp)

    return np.sum(np.square(rcnp-scnp))

def call_CN(Y):
    Y = Y.T
    Y_nor = []
    X = np.arange(min_ploidy, max_ploidy, 0.05)
    for cell in Y:
        RCNP = cell / np.mean(cell)
        sose_list = []
        for i in X:
            sose_list.append(sose(RCNP * i))
        ploidy = X[sose_list.index(min(sose_list))]
        SCNP = RCNP * ploidy
        Y_nor.append(SCNP)
    Y_nor = np.array(Y_nor)
    Y_nor = np.round(Y_nor)
    for i in range(len(Y_nor)):
        for j in range(len(Y_nor[0])):
            if Y_nor[i][j] > 10:
                Y_nor[i][j] = 10

    return Y_nor.T

def Gaussian_kernel(data,mu,sig):
    temp = np.array(-(data - mu) ** 2 / (2 * sig ** 2), dtype="float")
    density = np.exp(temp)

    return density

def get_norm_cell(cov_data):
    cov_std = stats.variation(cov_data, axis=0)
    kde = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(cov_std.reshape(-1, 1))
    dens = np.exp(kde.score_samples(np.arange(0, 1, 0.001).reshape(-1, 1)))
    peaks_pos = signal.argrelextrema(dens, np.greater)[0]
    first_peak = np.arange(0, 1, 0.001)[peaks_pos][0]
    density_list = Gaussian_kernel(cov_std, first_peak, 0.01)
    norm_index = []
    abnorm_index = []
    for i in range(len(density_list)):
        if density_list[i] > 1e-3:
            norm_index.append(i)
        else:
            abnorm_index.append(i)
    return cov_std, norm_index, abnorm_index

def Bias_norm(Y, norm_cell_index):
    Y =Y.T
    norm_cell_Y = Y[norm_cell_index]
    bias_matrix = []
    for cell in norm_cell_Y:
        bias_list = []
        median = np.median(cell)
        for bin in cell:
            bias = bin/median
            bias_list.append(bias)
        bias_list = np.array(bias_list)
        bias_matrix.append(bias_list)
    bias_matrix = np.array(bias_matrix)
    ave_bias = bias_matrix.mean(axis=0)
    ave_bias = np.where(ave_bias==0,1,ave_bias)
    gc_nor_Y = Y / ave_bias

    return gc_nor_Y.T

class fine_tune:
    def __init__(self, raw_cnv, Y):
        self.raw_cnv = raw_cnv.T
        self.Y = Y.T

    def cal_cn(self, alpha, read_num):
        p_list = []
        for copy_num in range(1, 8):
            p_list.append(stats.poisson.pmf(read_num, alpha * copy_num))
        check = 0
        for i in p_list:
            if i != 0:
                check = 1
                break
        if check == 0:
            cn = min(int(round(read_num / alpha)), 10)
        else:
            cn = p_list.index(max(p_list)) + 1

        return cn

    def cal_alpha(self, cn_list, read_list):
        alpha = np.sum(read_list) / np.sum(cn_list)
        return alpha

    def tune_cn(self):

        self.tuned_cn = []
        for i in range(len(self.Y)):
            reads = self.Y[i]
            cn = self.raw_cnv[i]
            diploid_index = []
            for j in range(len(cn)):
                if cn[j] == 2:
                    diploid_index.append(j)
            initial_alpha = np.mean(reads[diploid_index])/2
            alpha = initial_alpha
            reads_nor = []
            for i in reads:
                if i > alpha * 10:
                    reads_nor.append(alpha * 10)
                else:
                    reads_nor.append(i)
            reads_nor = np.array(reads_nor)

            while (1):
                cell_cn_list = []
                for read_num in reads_nor:
                    cell_cn = self.cal_cn(alpha, read_num)
                    cell_cn_list.append(cell_cn)

                old_alpha = alpha
                alpha = self.cal_alpha(cell_cn_list, reads_nor)

                if alpha == old_alpha:
                    break
            self.tuned_cn.append(np.array(cell_cn_list))
        return np.array(self.tuned_cn).T

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



class random_walk:
    def __init__(self, matrix):
        self.B_matrix = matrix
        self.M_1 = np.zeros((len(self.B_matrix), len(self.B_matrix[0])))
        for i in range(len(self.B_matrix)):
            sum_i = sum(self.B_matrix[i])
            if sum_i !=0:
                for j in range(len(self.B_matrix[0])):
                    self.M_1[i][j] = self.B_matrix[i][j] / sum_i
        self.M_1 = self.M_1.T

        self.M_2 = np.zeros((len(self.B_matrix), len(self.B_matrix[0])))
        for i in range(len(self.B_matrix[0])):
            sum_i = sum(self.B_matrix[:,i])
            if sum_i != 0:
                for j in range(len(self.B_matrix)):
                    self.M_2[j][i] = self.B_matrix[j][i] / sum_i

        self.transfer_M = np.dot(self.M_2, self.M_1)
        self.primary_d = np.diag(np.sum(self.B_matrix,axis=1))
        self.weight_matrix = np.dot(self.primary_d, self.transfer_M.T)
        self.p_list = np.diagonal(self.weight_matrix)

    def Gaussian_kernel(self, gamma):

        adj_matrix_1 = []
        for i in range(len(self.p_list)):
            dis_sq_matrix = (self.weight_matrix[i] - self.p_list[i]) ** 2
            adj_matrix_1.append(exp(-gamma * dis_sq_matrix))
        adj_matrix_1 = np.array(adj_matrix_1)

        adj_matrix_2 = []
        for i in range(len(self.p_list)):
            dis_sq_matrix = (self.weight_matrix[:, i] - self.p_list[i]) ** 2
            adj_matrix_2.append(exp(-gamma * dis_sq_matrix))
        adj_matrix_2 = np.array(adj_matrix_2)
        adj_matrix_2 = adj_matrix_2.T

        adj_matrix = (adj_matrix_1 + adj_matrix_2) / 2
        adj_matrix = KR_norm_juicer.KR_norm(adj_matrix)

        return adj_matrix

    def Gaussian_kernel_local(self, k):
        adj_matrix_1 = []
        for i in range(len(self.p_list)):
            dis_sq_matrix = (self.weight_matrix[i] - self.p_list[i]) ** 2
            dis_sq_matrix_filter = np.array(list(filter(lambda x: x != 0, dis_sq_matrix)))
            if len(dis_sq_matrix_filter) < k:
                k_dis = 1e-10
            else:
                k_dis = dis_sq_matrix_filter[np.argsort(dis_sq_matrix_filter)[k]]
            adj_matrix_1.append(exp(-dis_sq_matrix / k_dis))
        adj_matrix_1 = np.array(adj_matrix_1)
        adj_matrix_2 = []
        for i in range(len(self.p_list)):
            dis_sq_matrix = (self.weight_matrix[:, i] - self.p_list[i]) ** 2
            dis_sq_matrix_filter = np.array(list(filter(lambda x: x != 0, dis_sq_matrix)))
            if len(dis_sq_matrix_filter) < k:
                k_dis = 1e-10
            else:
                k_dis = dis_sq_matrix_filter[np.argsort(dis_sq_matrix_filter)[k]]
            adj_matrix_2.append(exp(-dis_sq_matrix / k_dis))
        adj_matrix_2 = np.array(adj_matrix_2)
        adj_matrix_2 = adj_matrix_2.T

        adj_matrix = (adj_matrix_1 + adj_matrix_2) / 2
        adj_matrix = KR_norm_juicer.KR_norm(adj_matrix)

        return adj_matrix

def average_seg(boundaries, chrom_matrix):
    chrom_matrix = chrom_matrix.T
    for i in range(len(chrom_matrix)):
        for j in range(len(boundaries)):
            if j == 0:
                chrom_matrix[i][0:boundaries[j]] = np.median(chrom_matrix[i][0:boundaries[j]])
            else:
                chrom_matrix[i][boundaries[j-1]:boundaries[j]] = np.median(chrom_matrix[i][boundaries[j-1]:boundaries[j]])

    return chrom_matrix.T


def save_matrix(matrix, bin_list, sample_list, var, cnv_name, meta_name):
    ploidy = np.average(matrix, axis=1)
    ploidy = np.round(ploidy, 1)
    df = pd.DataFrame(matrix, columns=bin_list, index=sample_list)
    df.to_csv(cnv_name)
    df_2 = pd.DataFrame(np.vstack((sample_list, ploidy, var)).T, columns=["cell_id", "c_ploidy", "c_variance"])
    df_2.to_csv(meta_name, index=False)

def main():

    Y, chr_name, bin_list, sample_list = read_matrix(os.path.join(work_path, "genome_cov.bed"))
    var, norm_cell_index, abnorm_cell_index = get_norm_cell(Y)
    cov_matrix = Bias_norm(Y, norm_cell_index)
    chosen_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    for chrom in chosen_chr:
        print("process %s..."%chrom)
        indexes = np.where(chr_name == chrom)
        chrom_matrix = cov_matrix[indexes]
        matrix = chrom_matrix[:, abnorm_cell_index]
        random_walk_sample = random_walk(matrix)
        adj_matrix = random_walk_sample.Gaussian_kernel_local(K)
        dp_process = DP_process(adj_matrix, 20)
        dp_process.dp_init()
        dp_process.dp_process()
        dp_process.find_k()
        boundaries = dp_process.find_boundaries()
        new_chrom = average_seg(boundaries, chrom_matrix)
        if chrom == "chr1":
            new_cov_matrix = new_chrom
        else:
            new_cov_matrix = np.vstack((new_cov_matrix, new_chrom))
    cn_matrix = call_CN(new_cov_matrix)
    save_matrix(cn_matrix.T, bin_list, sample_list, var, os.path.join(work_path, "cnv_matrix.csv"), os.path.join(work_path, "cnv_meta.csv"))


if __name__ == "__main__":
    main()

