import sys
import os
import linecache
import numpy as np

work_path = sys.argv[1]
ref = sys.argv[2]
bam_files = os.path.join(sys.argv[3], sys.argv[4])
bin_size = int(sys.argv[5])
ref_type = sys.argv[6]

chosen_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

def construct_bins(ref, bin_size):
    print("Constructing genome-wide consecutive bins...")
    genome_size = os.path.join(work_path, "genome_size.txt")
    genome_consecutive_bins = os.path.join(work_path, "genome_consecutive_bins.bed")
    genome_consecutive_bins_add = os.path.join(work_path, "genome_consecutive_bins_add.bed")
    os.system("faidx %s -i chromsizes > %s"%(ref, genome_size))
    f_in = open(genome_size)
    f_out = open(genome_consecutive_bins, "w")
    count = 1
    for line in f_in:
        chr_num = line.split("\t")[0]
        length = int(line.split("\t")[1])
        start = 1
        if chr_num in chosen_chr:
            while (length-start) > bin_size-1:
                f_out.write(chr_num)
                f_out.write("\t")
                f_out.write(str(start))
                f_out.write("\t")
                f_out.write(str(start+bin_size-1))
                f_out.write("\n")
                start = start + bin_size
                count += 1
            f_out.write(chr_num)
            f_out.write("\t")
            f_out.write(str(start))
            f_out.write("\t")
            f_out.write(str(length))
            f_out.write("\n")
            count += 1
    f_out.close()
    f_in.close()

    f_in = open(genome_size)
    f_out = open(genome_consecutive_bins_add, "w")
    count = 1
    for line in f_in:
        chr_num = line.split("\t")[0]
        length = int(line.split("\t")[1])
        start = 1
        if chr_num in chosen_chr:
            while (length-start) > bin_size-1:
                f_out.write(chr_num)
                f_out.write("\t")
                f_out.write(str(start))
                f_out.write("\t")
                f_out.write(str(start+bin_size-1))
                f_out.write("\t")
                f_out.write(str(count))
                f_out.write("\n")
                start = start + bin_size
                count += 1
            f_out.write(chr_num)
            f_out.write("\t")
            f_out.write(str(start))
            f_out.write("\t")
            f_out.write(str(length))
            f_out.write("\t")
            f_out.write(str(count))
            f_out.write("\n")
            count += 1
    f_out.close()
    f_in.close()

def cal_gc(bed_file):
    print("Calculating GC content...")
    genome_gc = os.path.join(work_path, "genome_gc.bed")
    os.system("bedtools nuc -fi %s -bed %s | cut -f 1-3,5 > %s"%(ref, bed_file, genome_gc))

def cal_map(bed_file):
    print("Calculating mappability...")
    genome_map = os.path.join(work_path, "genome_mappability.tab")
    if ref_type == "hg19":
        os.system("bigWigAverageOverBed hg19_mappability.bigWig %s %s"%(bed_file, genome_map))
    elif ref_type == "hg38":
        os.system("bigWigAverageOverBed hg38_mappability.bigWig %s %s"%(bed_file, genome_map))

def cal_cov(bed_file, bam_files):
    read_num = os.path.join(work_path, "read_num.txt")
    sample_list_f = os.path.join(work_path, "sample_list.txt")
    genome_cov_raw = os.path.join(work_path, "genome_cov_raw.bed")
    genome_cov = os.path.join(work_path, "genome_cov.bed")
    if os.path.exists(read_num):
        os.system("rm %s"%read_num)
    os.system("ls %s > %s"%(bam_files, sample_list_f))
    sample_list = []
    f = open(sample_list_f)
    for line in f:
        sample_list.append(line.strip("\n"))
    f.close()
    sample_list = np.array(sample_list)
    for sample in sample_list:
        os.system("samtools view -c %s >>%s"%(sample, read_num))
    read_num_list = []
    f = open(read_num)
    for line in f:
        read_num_list.append(int(line.strip("\n")))
    f.close()
    read_num_list = np.array(read_num_list)
    print("Calculating coverage...")
    os.system('/bin/bash -c "bedtools multicov -bams %s -bed %s -q 60 >%s"'%(bam_files, bed_file, genome_cov_raw))
    raw_cov = []    
    loc = []
    f = open(genome_cov_raw)
    for line in f:
        raw_cov.append(line.strip("\n").split("\t")[3:])
        loc.append(line.strip("\n").split("\t")[:3])
    f.close()
    raw_cov = np.array(raw_cov)
    loc = np.array(loc)
    raw_cov = raw_cov.astype(int)
    map_num_list = raw_cov.sum(axis=0)
    ratio_list = map_num_list/read_num_list
    remove_list = []
    for i in range(len(ratio_list)):
        if ratio_list[i] < 0.2:
            remove_list.append(i)
    remove_list.reverse()
    for index in remove_list:  
        sample_list = np.delete(sample_list,index,0)
        raw_cov = np.delete(raw_cov,index,1)
    f = open(genome_cov, "w")
    f.write("chromosome")
    f.write("\t")
    f.write("start")
    f.write("\t")
    f.write("stop")
    f.write("\t")
    for sample in sample_list:
        f.write(sample.split("/")[-1].split(".")[0])
        f.write("\t")
    f.write("\n")
    for i in range(raw_cov.shape[0]):
        f.write(loc[i,0])
        f.write("\t")
        f.write(loc[i,1])
        f.write("\t")
        f.write(loc[i,2])
        f.write("\t")
        for j in range(raw_cov.shape[1]):
            f.write(str(raw_cov[i,j]))
            f.write("\t")
        f.write("\n")
    f.close()


def filter_bins(gc_file, map_file):
    ref = os.path.join(work_path, "ref.bed")
    bin_num = len(open(map_file, 'r').readlines())
    f_out = open(ref, "w")
    count = 0
    for i in range(bin_num):
        gc_line = linecache.getline(gc_file, i+2).strip("\n").split("\t")
        gc = float(gc_line[3])
        mappability = float(linecache.getline(map_file, i+1).strip("\n").split("\t")[5])
        if gc > 0.2 and gc < 0.8 and mappability > 0.9:
            count += 1
            f_out.write(gc_line[0])
            f_out.write("\t")
            f_out.write(gc_line[1])
            f_out.write("\t")
            f_out.write(gc_line[2])
            f_out.write("\n")
    f_out.close()

def main():
    construct_bins(ref, bin_size)
    cal_gc(os.path.join(work_path, "genome_consecutive_bins.bed"))
    cal_map(os.path.join(work_path, "genome_consecutive_bins_add.bed"))
    filter_bins(os.path.join(work_path, "genome_gc.bed"), os.path.join(work_path, "genome_mappability.tab"))
    cal_cov(os.path.join(work_path, "ref.bed"), bam_files)

if __name__ == "__main__":
    main()
