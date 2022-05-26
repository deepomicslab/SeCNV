"""SeCNV
Usage:
    SeCNV.py <bam_path> <output_path> <ref_file> [--ref=<rn>] [--bin_size=<bl>] [--min_ploidy=<p>] [--max_ploidy=<p>] [--pattern=<pn>] [--topK=<n>] [--sigma=<n>] [--normal_cell=<n>]
    SeCNV.py (-h | --help)
    SeCNV.py --version

Options:
    -r --ref=<rn>   The reference used (hg19 or hg38) [default: hg19].
    -b --bin_size=<bl>  The length of bin [default: 500000].
    -min --min_ploidy=<p>  The minimal ploidy [default: 1.5].
    -max --max_ploidy=<p>  The maximal ploidy [default: 5.0].
    -p --pattern=<pn>    The pattern of bam file names [default: *dedup.bam].
    -K --topK=<n>  The K largest distances used to construct adjacency matrix [default: auto_set].
    -s --sigma=<n>  The standard deviation of the Gaussian kernel function [default: auto_set].
    -n --normal_cell=<n>  The file with normal cell IDs [default: None].
    -h --help   Show this screen.
    -v --version    Show version.
"""
from docopt import docopt
import os
import time

def main(arguments):
    start_time = time.time()
    bam_path = arguments.get("<bam_path>")
    output_path = arguments.get("<output_path>")
    ref = arguments.get("<ref_file>")
    ref_type = arguments.get("--ref")
    bin_size = arguments.get("--bin_size")
    min_ploidy = arguments.get("--min_ploidy")
    max_ploidy = arguments.get("--max_ploidy")
    bam_pattern = arguments.get("--pattern")
    K = str(arguments.get("--topK"))
    sigma = str(arguments.get("--sigma"))
    normal_cell = str(arguments.get("--normal_cell"))
    os.system("python preprocess.py %s %s %s %s %s %s"%(output_path, ref, bam_path, bam_pattern, bin_size, ref_type))
    os.system("python call_cn.py %s %s %s %s %s %s"%(output_path, min_ploidy, max_ploidy, K, sigma, normal_cell))
    stop_time = time.time()
    print("Finish. Cost time %s s"%(stop_time-start_time))

if __name__ == "__main__":
    arguments = docopt(__doc__, version="SeCNV 0.1.1")
    main(arguments)
