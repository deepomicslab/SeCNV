"""SeCNV
Usage:
    SeCNV.py <bam_path> <output_path> <ref_file> [--ref=<rn>] [--bin_size=<bl>] [--min_ploidy=<p>] [--max_ploidy=<p>] [--pattern=<pn>] [--K_neighbor=<n>]
    SeCNV.py (-h | --help)
    SeCNV.py --version

Options:
    -r --ref=<rn>   The reference used (hg19 or hg38) [default: hg19].
    -b --bin_size=<bl>  The length of bin [default: 500000].
    -min --min_ploidy=<p>  The length of bin [default: 1.5].
    -max --max_ploidy=<p>  The length of bin [default: 5.0].
    -p --pattern=<pn>    The pattern of bam file names [default: *dedup.bam].
    -K --K_neighbor=<n>	The number of neighbors [default: 5].
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
    K = arguments.get("--K_neighbor")
    os.system("python preprocess.py %s %s %s %s %s %s"%(output_path, ref, bam_path, bam_pattern, bin_size, ref_type))
    os.system("python call_cn.py %s %s %s %s"%(output_path, min_ploidy, max_ploidy, K))
    stop_time = time.time()
    print(stop_time-start_time)

if __name__ == "__main__":
    arguments = docopt(__doc__, version="SeCNV 0.1.0")
    main(arguments)
