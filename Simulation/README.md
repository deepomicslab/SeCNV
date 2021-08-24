# Simulation experiment
+ simulate\_data.py: The script to simulate read coverge matrix with different breakpoint distributions and varying amplitudes of noises. Usage: python simulate\_data.py norm\_cell\_percentage breakpoint\_occurrence\_probability gaussian\_noise\_variance. Output: data\_norm\_bp\_gn.txt:read coverage matrix (bin\_num * cell\_num);data\_norm\_bp\_gn.png:copy number profiles visualization.

+ SeCNV\_seg.py: The segmentation script for SeCNV. Usage: python SeCNV\_seg.py data\_norm\_bp\_gn.txt
