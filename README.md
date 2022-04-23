# SeCNV
SeCNV, a single cell copy number profiling tool.

![image](https://github.com/deepomicslab/SeCNV/blob/main/Framework.jpg)

## Prerequisite
Please install Git LFS and use the following command to clone the repository.
```shell
git lfs clone https://github.com/deepomicslab/SeCNV.git
```

The scripts are written in Python3. Following Python packages should be installed:
+ numpy
+ pandas
+ scipy
+ sklearn
+ linecache

Since there are some Bioinformatic pre-processing steps in SeCNV pipeline, following Bioinformatic tools shoud be installed and set environment variables.
+ bwa
+ samtools
+ bedtools
+ bigWigAverageOverBed
+ pyfaidx
+ picard 

## Data preparation
To run SeCNV, the bigwig files, hg19\_mappability.bigWig and hg38\_mappability.bigWig, should be downloaded from [Google Drive](https://drive.google.com/drive/folders/1XGuXa9muRtOiAjtz1J4gWO5Qk3c5-PwS?usp=sharing) and put under the Script folder.

The reference hg19 or hg38, which can be downloaded from NCBI, should be prepared and built index. 

## Usage
### Bioinformatic pre-processing
For alignment, sorting, adding read group, and deduplication, we recommend the following steps.
+ Align FASTQ files to the reference genome
```shell
bwa mem -M -t 8 hg19.fa file_name.fastq.gz > file_name.sam
samtools view -bS file_name.sam >file_name.bam 
```
+ Sort
```shell
java -Xmx30G -jar picard.jar SortSam INPUT=file_name.bam OUTPUT=file_name.sorted.bam SORT_ORDER=coordinate 
```
+ Add read Group
```shell
java -Xmx40G -jar picard.jar AddOrReplaceReadGroups I=file_name.sorted.bam O=file_name.sorted.rg.bam RGID=file_name RGLB=NAVIN_Et_Al RGPL=ILLUMINA RGPU=machine RGSM=file_name
```
+ Dedup
```shell
java -Xmx40G -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=file_name.sorted.rg.bam O=file_name.sorted.rg.dedup.bam METRICS_FILE=file_name.sorted.rg.dedup.metrics.txt PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_VERSION=null PROGRAM_GROUP_NAME=MarkDuplicates
java -jar picard.jar BuildBamIndex I=file_name.sorted.rg.dedup.bam
```
Please change **hg19.fa** to your reference location and **file\_name** to your FASTQ file name.
 
### SeCNV
Next, SeCNV takes the bam files as input to profile copy number.
```shell
cd Scripts
python SeCNV.py input_fold output_fold ref_file
```
Input\_fold is where the bam files are, output\_file is where the output files will be (an empty fold is recommended), and the ref\_file is the path of the indexed reference hg19 or hg38. Other parameters are shown bellow:

+ -r or --ref:	The reference used (hg19 or hg38) [default: hg19].
+ -b or --bin\_size:	The length of bin [default: 500000].
+ -min or --min\_ploidy:	The minimal ploidy [default: 1.5].
+ -max or --max\_ploidy:	The maximal ploidy [default: 5].
+ -p or --pattern:    The pattern of bam file names [default: \*dedup.bam]. 
+ -K or --topK:	The K largest distances used to construct adjacency matrix [default: auto\_set].
+ -s or --sigma:	The standard deviation of the Gaussian kernel function [default: auto\_set].

For more information, please use python SeCNV.py -h or python SeCNV.py --help.

### Maintainer
WANG Ruohan ruohawang2-c@my.cityu.edu.hk
