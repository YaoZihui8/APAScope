import argparse

parser = argparse.ArgumentParser(description="RNA-seq APA analysis pipeline.")
parser.add_argument("-f", "--fastq", required=True,
                    help="input raw RNA-seq data directory (raw data in fastq format) e.g. /public5/lisj/00_raw. Attention: Please name the raw data as *_1.fastq and *_2.fastq")
# parser.add_argument("-r", "--reference", type=str,
#                     help="reference genome fasta file e.g. /public5/lisj/Genome/fasta.fna")
parser.add_argument("-o", "--output", required=True, help="output directory e.g. /public5/lisj/ ")
parser.add_argument("-rl", "--readlimit", type=str, required=False,default = '20', help=" ")
parser.add_argument("-eg", "--endgap", type=str, required=False,default = '10', help=" ")
# parser.add_argument("-c", "--chr", required=True, help=""
# parser.add_argument("-s", "--sex", required=True, help="")

args = parser.parse_args()

print("checking modules")
import glob, os, sys, time
import pandas as pd
import numpy as np
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects import StrVector
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects

print("imports done")

# reference_dir = args.reference.rstrip('/')
library_dir = os.path.dirname(os.path.abspath(__file__))
tools_dir = os.path.join(library_dir, 'tools')

all_result_dir = os.path.join(args.output, 'output')
os.makedirs(all_result_dir, exist_ok=True)

### 01 fastp ###
in_dir = args.fastq.rstrip('/')
out_dir = os.path.join(all_result_dir, '01_CleanData')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)
for filename in glob.glob(in_dir + '/*_1.fastq.gz'):
    basename = os.path.basename(filename)
    id = basename.split('_1.fastq.gz')[0]
    if os.path.isfile(id + '_fastp.html') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' Trimming: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')
        # cmd = (
        #         'fastp -w 6 -i ' + filename + ' -I ' + in_dir + '/' + id + '_2.fastq.gz -o ' + id + '_1.fastq.gz -O ' + id
        #         + '_2.fastq.gz -f 10 -t 2 -q 20 -u 20 -e 20 -n 0 -l 75 -y -c -j ' + id + '_fastp.json -h ' + id + '_fastp.html -R "' + id + ' fastp report" 2> ' + id + '_fastp.log')
        cmd = (
                'fastp -w 16 -i ' + filename + ' -I ' + in_dir + '/' + id + '_2.fastq.gz -o ' + id + '_1.fastq.gz -O ' + id
                + '_2.fastq.gz -h ' + id + '_fastp.html -R "' + id + ' fastp report" 2> ' + id + '_fastp.log')
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('FASTQ already trimmed before: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')
os.chdir(all_result_dir)

## HISAT2 INDEX
# if os.path.isfile(reference_dir + '/genomic.8.ht2') == False:
#     os.chdir(reference_dir)
#     cmd = 'hisat2-build -p 8 ' + reference_dir + '/genomic.fa ' + 'genomic'
#     print('Building HISAT2 Index!')
#     print(cmd)
#     os.system(cmd)
# else:
#     print('HISAT2 index has already been built')
# os.chdir(all_result_dir)

### Align to genome ###
in_dir = os.path.join(all_result_dir, '01_CleanData')
out_dir = os.path.join(all_result_dir, '02_HISAT2')
out_dir2 = os.path.join(all_result_dir, 'B01_Unmapped')
os.makedirs(out_dir, exist_ok=True)
os.makedirs(out_dir2, exist_ok=True)
os.chdir(out_dir)
for filename in glob.glob(in_dir + '/*_1.fastq.gz'):
    basename = os.path.basename(filename)
    id = basename.split('_1.fastq.gz')[0]
    if os.path.isfile(out_dir + '/' + id + '.align.log') == False:
        print(time.asctime(time.localtime(
            time.time())) + ' Align to genome by HISAT2: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')
        cmd = 'hisat2 -p 40 --dta -x ' + library_dir + '/hg19/hg19 -1 ' + in_dir + '/' + id + '_1.fastq.gz -2 ' + in_dir + '/' + id + '_2.fastq.gz --rna-strandness RF --no-mixed --no-discordant --un-conc-gz ' + out_dir2 + '/' + id + '_um_%.fastq.gz --new-summary --summary-file ' + out_dir + '/' + id + '.align.log --no-unal | samtools sort -@ 9 -m 5G -o ' + out_dir + '/' + id + '.bam'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('Aligning already done before: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')
os.chdir(all_result_dir)

### switch to fasta ###
in_dir = os.path.join(all_result_dir, '02_HISAT2')
out_dir = os.path.join(all_result_dir, '03_fastq')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)
for filename in glob.glob(in_dir + '/*.bam'):
    basename = os.path.basename(filename)
    id = os.path.splitext(basename)[0]
    if os.path.isfile(out_dir + '/' + id + '.fa') == False:
        print(time.asctime(time.localtime(
            time.time())) + ' Switch to fasta: ' + id + '.bam')
        cmd = 'samtools fasta ' + in_dir + '/' + id + '.bam > ' + id + '.fa --threads 40'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('Fasta already done before: ' + id + '.bam')
os.chdir(all_result_dir)

### 04 obtain PAS database ###
PAS_dir = os.path.join(all_result_dir, '04_PAS_database')
os.makedirs(PAS_dir, exist_ok=True)

######## 01 FindTail ###
in_dir = os.path.join(all_result_dir, '03_fastq')
out_dir = os.path.join(PAS_dir, '01_FindTail')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.fa'):
    basename = os.path.basename(filename)
    id = os.path.splitext(basename)[0]
    if os.path.isfile(id + '.pr.pA') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' FindTail: ' + id + '.fa')
        cmd = f'{tools_dir}/findtail_v1.01 --input_file ' + in_dir + '/' + id + f'.fa --seqlength 1000 --endgap {args.endgap} --taillength 15 --identity 85 --ptype A --stype T --output_format fasta --output_type pr --d > ' + out_dir + '/' + id + '.pr.pA'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('FindTail already done before: ' + id + '.pr.pA')
os.chdir(all_result_dir)

######## 02 FilterIdentity ###
in_dir = os.path.join(PAS_dir, '01_FindTail')
out_dir = os.path.join(PAS_dir, '02_FilterIdentity')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.pr.pA'):
    basename = os.path.basename(filename)
    id = basename.split('.pr.pA')[0]
    if os.path.isfile(id + '.i94') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' FilterIdentity: ' + id + '.pr.pA')
        cmd = tools_dir + '/FilterIdentity --input_file ' + in_dir + '/' + id + '.pr.pA --seqlength 500 --identity 94 --taillength 15 > ' + out_dir + '/' + id + '.i94'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('FilterIdentity already done before: ' + id + '.i94')
os.chdir(all_result_dir)

######## 03 ReplaceN ###
in_dir = os.path.join(PAS_dir, '02_FilterIdentity')
out_dir = os.path.join(PAS_dir, '03_ReplaceN')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.i94'):
    basename = os.path.basename(filename)
    id = basename.split('.i94')[0]
    if os.path.isfile(id + '.RN.i94') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' ReplaceN: ' + id + '.pr.pA')
        cmd = tools_dir + '/ReplaceN --input_file ' + in_dir + '/' + id + '.i94 --seqlength 500 --identity 94 --Ncount 5 > ' + out_dir + '/' + id + '.RN.i94'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('ReplaceN already done before: ' + id + '.RN.i94')
os.chdir(all_result_dir)

######## 04 TrimPolyA ###
in_dir = os.path.join(PAS_dir, '03_ReplaceN')
out_dir = os.path.join(PAS_dir, '04_TrimPolyA')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.RN.i94'):
    basename = os.path.basename(filename)
    id = basename.split('.RN.i94')[0]
    if os.path.isfile(id + '.RN.fasta') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' TrimPolyA: ' + id + '.RN.i94')
        cmd = tools_dir + '/TrimPolyA --input_file ' + in_dir + '/' + id + '.RN.i94 --seqlength 500 --identity 94 --trimside r > ' + out_dir + '/' + id + '.RN.fasta'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('TrimPolyA already done before: ' + id + '.RN.fasta')
os.chdir(all_result_dir)

# ## Bowtie2 INDEX
# if os.path.isfile(reference_dir + '/genomic.4.bt2') == False:
#     os.chdir(reference_dir)
#     cmd = 'bowtie2-build -p 8 ' + reference_dir + '/genomic.fa ' + 'genomic'
#     print('Building Bowtie2 Index!')
#     print(cmd)
#     os.system(cmd)
# else:
#     print('Bowtie2 index has already been built')
# os.chdir(all_result_dir)

######## 05 Bowtie2 ###
in_dir = os.path.join(PAS_dir, '04_TrimPolyA')
out_dir = os.path.join(PAS_dir, '05_Bowtie2')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.RN.fasta'):
    basename = os.path.basename(filename)
    id = basename.split('.RN.fasta')[0]
    if os.path.isfile(id + '.sam') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' Bowtie2: ' + id + '.RN.fasta')
        cmd = 'bowtie2 -p 32 -f -x ' + library_dir + '/hg19/hg19 -U ' + in_dir + '/' + id + '.RN.fasta --no-unal -S ' + out_dir + '/' + id + '.sam'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('Bowtie2 already done before: ' + id + '.sam')
os.chdir(all_result_dir)

######## 06 CleanSam ###
in_dir = os.path.join(PAS_dir, '05_Bowtie2')
out_dir = os.path.join(PAS_dir, '06_CleanSam')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

for filename in glob.glob(in_dir + '/*.sam'):
    basename = os.path.basename(filename)
    id = basename.split('.sam')[0]
    if os.path.isfile(id + '.sam') == False:
        print(time.asctime(
            time.localtime(time.time())) + ' CleanSam: ' + id + '.sam')
        cmd = "grep 'AS:' " + in_dir + '/' + id + ".sam | grep -v 'XS:' > " + out_dir + '/' + id + '.sam'
        print('CMD: ' + cmd)
        os.system(cmd)
        cmd = f"perl {tools_dir}/cleansam_pr_pa_hg19.pl " + id + '.sam'
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print('CleanSam already done before: ' + id + '.sam')

######## 07 FilterPolyASite ###
in_dir = os.path.join(PAS_dir, '06_CleanSam')
out_dir = os.path.join(PAS_dir, '07_FilterPolyASite')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)
cmd = "cat " + in_dir + '/' + "*clean > " + out_dir + '/' + 'pA.sam.clean'
print('CMD: ' + cmd)
os.system(cmd)

# if args.chr == 1:
#     chromosomes = [str(x) for x in range(1, int(args.chr))] + ['X', 'Y']
# else:
#     chromosomes = [str(x) for x in range(1, int(args.chr))]
if os.path.isfile('highqc_polya_cluster.txt') == False:
    chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']
    for chrom in chromosomes:
        # 正链
        cmd = tools_dir + f"/FilterGenome --input_file pA.sam.clean --seqlength 500 --strand 0 --chromosome_name chr{chrom} > chr{chrom}_positive.fasta"
        print(f"CMD: {cmd}")
        os.system(cmd)

        # 负链
        cmd = tools_dir + f"/FilterGenome --input_file pA.sam.clean --seqlength 500 --strand 16 --chromosome_name chr{chrom} > chr{chrom}_negative.fasta"
        print(f"CMD: {cmd}")
        os.system(cmd)
    for fa_file in glob.glob("chr*.fasta"):
        output_file = fa_file.replace('.fasta', '_3reads.fasta')
        cmd = (
            f"{tools_dir}/FilterPolyASite "
            f"--input_file {fa_file} --seqlength 500 --read_limit 3 > {output_file}"
        )
        print(f"CMD: {cmd}")
        os.system(cmd)
    for chrom in chromosomes:
        for strand in ['positive', 'negative']:
            input_file = f"chr{chrom}_{strand}.fasta"
            output_file = input_file.replace('.fasta', '_cluster.fasta')
            cmd = (
                f"{tools_dir}/GetPolyaSiteCluster "
                f"--input_file {input_file} --seqlength 500 --read_limit {args.readlimit} "
                f"--cluster_seed 24 > {output_file}"
            )
            print(f"CMD: {cmd}")
            os.system(cmd)

    cmd = "cat *_cluster.fasta | sort | uniq > highqc_polya_cluster.txt"
    print(f"CMD: {cmd}")
    os.system(cmd)
else:
    print('highqc_polya_cluster.txt already done before')