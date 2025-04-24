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

### 05 QAPA ###
in_dir = os.path.join(PAS_dir, '07_FilterPolyASite')
out_dir = os.path.join(all_result_dir, '05_QAPA')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

cmd = r"""awk '{print $1"\t"$6-1"\t"$6"\t"$1":"$6":"$4":"$10"\t"$10"\t"$4}' """ + in_dir + r"""/highqc_polya_cluster.txt > pA.bed"""
print('CMD:', cmd)
os.system(cmd)

cmd = (
    f'conda run -n qapa qapa build --db {tools_dir}/qapa/examples/hg19/ensembl.identifiers.txt '
    f'-o pA.bed {library_dir}/hg19/gtf/gencode.v45lift37.annotation.genePred '
    '> pA_utrs.bed'
)
print('CMD:', cmd)
os.system(cmd)

cmd = r"""awk -F '\t' '{print $4}' pA_utrs.bed | awk -F '_' '{print $1}' | sed 's/"//g' > pA_utrs.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"paste -d '\t' utrs.bed pA_utrs.txt > pA_utrs_new.bed"
print('CMD:', cmd)
os.system(cmd)

# qapa_merge_bed.R
# 读取文件 (对应R的read.table)
file1 = pd.read_csv("pA_utrs_new.bed", sep='\t', header=None)
file2 = pd.read_csv(f"{library_dir}/qapa/examples/hg19/ensembl.identifiers.txt", sep='\t', header=None)

# 重命名列 (对应colnames)
file1 = file1.rename(columns={7: 'enst'})  # 第8列索引为7
file2 = file2.rename(columns={1: 'enst'})  # 第2列索引为1

# 合并数据 (对应merge)
merged_data = pd.merge(
    file1,
    file2,
    on='enst',
    how='left',  # all.x=TRUE
    suffixes=('.x', '.y')
)

# 选择列并重命名 (对应列选择)
data = merged_data[[
    '0.x', '1', '2.x', '3.x', '4.x', '5', '6', '4.y'
]].rename(columns={
    '0.x': 'V1.x',
    '1': 'V2',
    '2.x': 'V3.x',
    '3.x': 'V4.x',
    '4.x': 'V5.x',
    '5': 'V6',
    '6': 'V7',
    '4.y': 'V5.y'
})

# 条件替换逻辑 (对应mutate + if_else)
# 处理NaN值并转换类型
data['V4.x'] = data['V4.x'].astype(str)
data['V7'] = data['V7'].fillna('').astype(str)
data['V5.y'] = data['V5.y'].fillna('').astype(str)

# 使用正则表达式进行替换
mask = data.apply(lambda x: x['V7'] in x['V4.x'], axis=1)
data.loc[mask, 'V4.x'] = data[mask].apply(
    lambda x: x['V4.x'].replace(x['V7'], x['V5.y'], regex=False),
    axis=1
)

# 最终列选择 (对应最后的列筛选)
output_data = data[['V1.x', 'V2', 'V3.x', 'V4.x', 'V5.x', 'V6']]

# 输出结果 (对应write.table)
output_data.to_csv(
    "output.bed",
    sep='\t',
    header=False,
    index=False,
    quoting=3  # 相当于quote=FALSE
)

cmd = "sed 's/hsa/hg19/g' output.bed > pA_utrs.bed"
print('CMD:', cmd)
os.system(cmd)


def read_file(filename):
    data = set()
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if "\t" in line:
                chrom, start, end, utr, length, strand = line.split("\t")
                start = int(start)
                end = int(end)
                length = int(length)
                data.add((chrom, start, end, utr, length, strand))
    return data


new = read_file('pA_utrs.bed')
ref = read_file(f'{tools_dir}/qapa/examples/hg19/qapa_3utrs.gencode.hg19.bed')


def find_overlap(set1, set2):
    overlap = set()
    non_overlap = set()

    for item1 in set1:
        chrom1, start1, end1, utr1, length1, strand1 = item1
        found_overlap = False

        for item2 in set2:
            chrom2, start2, end2, utr2, length2, strand2 = item2

            if chrom1 == chrom2 and strand1 == strand2:
                if strand1 == '+':
                    if abs(start1 - start2) <= 25 and abs(length1 - length2) <= 25:
                        overlap.add(item1)
                        found_overlap = True
                        break
                elif strand1 == '-':
                    if abs(end1 - end2) <= 25 and abs(length1 - length2) <= 25:
                        overlap.add(item1)
                        found_overlap = True
                        break

        if not found_overlap:
            non_overlap.add(item1)

    return overlap, non_overlap


overlap, non_overlap = find_overlap(new, ref)


def save_overlap_to_bed(overlap, output_filename):
    with open(output_filename, "w") as file:
        for item in overlap:
            file.write("\t".join(map(str, item)) + "\n")


def save_non_overlap_to_bed(non_overlap, output_filename):
    with open(output_filename, "w") as file:
        for item in non_overlap:
            file.write("\t".join(map(str, item)) + "\n")


output_overlap_file = 'overlapping_utrs.bed'
output_non_overlap_file = 'non_overlapping_utrs.bed'

save_overlap_to_bed(overlap, output_overlap_file)
save_non_overlap_to_bed(non_overlap, output_non_overlap_file)

cmd = f"cat {tools_dir}/qapa/examples/hg19/qapa_3utrs.gencode.hg19.bed non_overlapping_utrs.bed | sort -k1,1 -k2,2n > qapa_3utrs_updata.gencode.hg19.bed"
print('CMD:', cmd)
os.system(cmd)

cmd = "grep '+' qapa_3utrs_updata.gencode.hg19.bed > P.bed"
print('CMD:', cmd)
os.system(cmd)

cmd = "grep -v '+' qapa_3utrs_updata.gencode.hg19.bed > N.bed"
print('CMD:', cmd)
os.system(cmd)

cmd = "cat P.bed N.bed > qapa_3utrs_updata.gencode.hg19.bed"
print('CMD:', cmd)
os.system(cmd)

cmd = f"conda run -n qapa qapa fasta -f {tools_dir}/qapa/examples/hg19/hg19.fa qapa_3utrs_updata.gencode.hg19.bed qapa_3utrs_updata.gencode.hg19.fa" + tools_dir + "/salmon-latest_linux_x86_64/bin/salmon index -t qapa_3utrs_updata.gencode.hg19.fa -i update_utr_library"
print('CMD:', cmd)
os.system(cmd)

### QAPA - Salmon###
in_dir2 = os.path.join(all_result_dir, '01_CleanData')
out_dir2 = os.path.join(out_dir, 'Salmon')
os.makedirs(out_dir2, exist_ok=True)
for filename in glob.glob(in_dir2 + '/*_1.fastq.gz'):
    basename = os.path.basename(filename)
    id = basename.split('_1.fastq.gz')[0]
    if os.path.isfile('align_and_estimate_abundance.log') == False:
        print(time.asctime(time.localtime(
            time.time())) + ' Run Salmon: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')
        cmd = (
            f"{tools_dir}/salmon-latest_linux_x86_64/bin/salmon quant -l A "
            f"-1 {in_dir2}/{id}_1.fastq.gz -2 {in_dir2}/{id}_2.fastq.gz "
            f"-i {tools_dir}/qapa/utr_library/ "
            f"-o Salmon/{id}.salmon.count "
            f"-p 20"
        )
        print('CMD: ' + cmd)
        os.system(cmd)
    else:
        print(' Run Salmon: ' + id + '_1.fastq.gz and ' + id + '_2.fastq.gz')

cmd = f"conda run -n qapa qapa quant --db {tools_dir}/qapa/examples/hg19/ensembl.identifiers.txt Salmon/*.salmon.count/quant.sf > pau_results.txt"
print('CMD:', cmd)
os.system(cmd)

### 06 APAlyzer ###
in_dir = os.path.join(PAS_dir, '07_FilterPolyASite')
out_dir = os.path.join(all_result_dir, '06_APAlyzer')
os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

cmd = r"""awk '{print $1"\t"$2"\t"$3"\t"$1":"$6":"$4"\t"$10"\t"$4}' """ + in_dir + """highqc_polya_cluster.txt > cancer28_PAS.bed"""
print('CMD:', cmd)
os.system(cmd)

cmd = f"bedtools intersect -a PAS.bed -b {library_dir}/hg19/gtf/gtf_transcript.bed -wa -wb -s > intersection_file_transcript.txt"
print('CMD:', cmd)
os.system(cmd)
cmd = r"awk 'BEGIN{FS=\" \";OFS=\"\t\"}{print $1,$2,$3,$15,$6,$10,$13,$14,$2}' intersection_file_transcript.txt > dfPAS.txt"
print('CMD:', cmd)
os.system(cmd)

numpy2ri.activate()
with localconverter(default_converter + numpy2ri.converter):
    robjects.r.source(library_dir+'/apalyzer.R')
    result = \
        robjects.r.apalyzer(out_dir)[0]


def read_file(filename):
    data = set()
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if ":" in line:
                PASid, gene, chrom, strand, pos, upSS, downSS = line.split("\t")
                pos = int(pos)
                data.add((PASid, gene, chrom, strand, pos, upSS, downSS))
    return data


REF = read_file(f'{tools_dir}/apalyzer/refIPA_hg19.tsv')
NEW = read_file('IPAraw.tsv')


def find_overlap(set1, set2):
    overlap = set()
    for item1 in set1:
        PASid1, gene1, chrom1, strand1, pos1, upSS1, downSS1 = item1
        for item2 in set2:
            PASid2, gene2, chrom2, strand2, pos2, upSS2, downSS2 = item2
            if (chrom1 == chrom2 and
                    strand1 == strand2 and
                    abs(pos1 - pos2) <= 25):
                overlap.add(item1)
                break
    return overlap


overlap = find_overlap(NEW, REF)

non_overlap_NEW = set(NEW) - overlap

with open('exist_pA.tsv', 'w') as overlap_file:
    overlap_file.write("PASid\tgene_symbol\tChrom\tStrand\tPos\tupstreamSS\tdownstreamSS\n")  # Header
    for item in overlap:
        overlap_file.write("\t".join(map(str, item)) + "\n")

with open('non_exist_pA.tsv', 'w') as non_overlap_file:
    non_overlap_file.write("PASid\tgene_symbol\tChrom\tStrand\tPos\tupstreamSS\tdownstreamSS\n")  # Header
    for item in non_overlap_NEW:
        non_overlap_file.write("\t".join(map(str, item)) + "\n")

print("Results have been written to exist_pA.tsv and non_exist_pA.tsv")

cmd = """grep 'ENSG' non_exist_pA.tsv | awk -F '\t' '{print $2}' | sed 's/"//g' | sort | uniq  > non_exist_pA_ensg.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = f"for i in $(cat non_exist_pA_ensg.txt); do echo -n $i"  f"; grep $i {library_dir}/hg19/gencode.v23.annotation.gene.probemap; echo ""; done | awk '{print $1,$3}' > non_exist_pA_ensg2symbol.txt"
print('CMD:', cmd)
os.system(cmd)

cmd = "sed -i '/^\s*$/d' non_exist_pA_ensg2symbol.txt"
print('CMD:', cmd)
os.system(cmd)

cmd = ""
print('CMD:', cmd)
os.system(cmd)

cmd = """while IFS= read -r line; do i=$(echo "$line" | awk '{print $1}') j=$(echo "$line" | awk '{print $2}') sed -i "s/${i}/${j}/" tcgaIPAraw.tsv done < non_exist_pA_ensg2symbol.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = """while IFS= read -r line; do
    i=$(echo "$line" | awk '{print $1}')
    j=$(echo "$line" | awk '{print $2}')
    sed -i "s/${i}/${j}/" non_exist_pA.tsv
done < non_exist_pA_ensg2symbol.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = "cat dfIPA_hg19.tsv updIPA/non_exist_pA.tsv | sort -k2,2 -k5,5n > updIPA_hg19.tsv"
print('CMD:', cmd)
os.system(cmd)

cmd = r"awk -F '\t' '{print $1}' dfLE_hg19.tsv | sort | uniq > updIPA/dfLE_gene.txt"
print('CMD:', cmd)
os.system(cmd)

cmd = r"awk -F '\t' '{print $2}' non_exist_pA.tsv | sort | uniq > tcgaLE_gene.txt"
print('CMD:', cmd)
os.system(cmd)

cmd = "for i in $(comm -23 LE_gene.txt dfLE_gene.txt); do grep $i tcgaLEraw.tsv; done > newLEraw.tsv"
print('CMD:', cmd)
os.system(cmd)

cmd = "cat dfLE_hg19.tsv updIPA/newLEraw.tsv | sort -k1,1n > updLE_hg19.tsv"
print('CMD:', cmd)
os.system(cmd)

cmd = r"""grep '+' UTRraw.tsv | awk '{print $2"\t"$6"\t"$5"\t"$4"\t"$1"\t"$3}' | sed 's/"//g' > UTRraw.bed"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"""grep '-' UTRraw.tsv | awk '{print $2"\t"$5"\t"$6"\t"$4"\t"$1"\t"$3}' | sed 's/"//g' >> UTRraw.bed"""
print('CMD:', cmd)
os.system(cmd)

cmd = f"bedtools intersect -a UTRraw.bed -b {library_dir}/apalyzer/refUTRraw_hg19.bed -v > UTR_new.bed"
print('CMD:', cmd)
os.system(cmd)

cmd = r"""grep '+' UTR_new.bed | awk -F '\t' '{print $5"\t"$1"\t"$6"\t"$4"\t"$3"\t"$2}' > newUTR.tsv"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"""grep -v '+' UTR_new.bed | awk -F '\t' '{print $5"\t"$1"\t"$6"\t"$4"\t"$2"\t"$3}' >> newUTR.tsv"""
print('CMD:', cmd)
os.system(cmd)

cmd = """awk '{if($4 != $5) print $0}' newUTR.tsv > processed_newUTR.tsv"""
print('CMD:', cmd)
os.system(cmd)

cmd = f"bedtools intersect -a {tools_dir}/apalyzer/refUTRraw_hg19.bed -b PAS.bed -s -wa -wb > update_pas.txt"
print('CMD:', cmd)
os.system(cmd)

cmd = """awk '{print $5}' update_pas.txt | sort | uniq -c | awk -v id='1' '$1 == id {print $2}' > update_1.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = """for i in $(cat update_1.txt); do awk -v gene="$i" '$5 == gene {print $0}' update_pas.txt; done | sed 's/:/\t/g' > update_1pas.txt"""
print('CMD:', cmd)
os.system(cmd)

input_file = 'update_1pas.txt'
overlap_lines = []
non_overlap_lines = []

with open(input_file, 'r') as infile:
    for line in infile:
        columns = line.strip().split('\t')
        if len(columns) >= 11:
            try:
                pA = int(columns[10])
                strand = columns[5]
                if strand == '+':
                    last = int(columns[2]) - 2000
                    first = int(columns[3])
                elif strand == '-':
                    last = int(columns[1]) + 2000
                    first = int(columns[3])
                else:
                    continue

                if abs(pA - first) <= 25 or abs(pA - last) <= 25:
                    overlap_lines.append(line)
                else:
                    non_overlap_lines.append(line)

            except ValueError:
                print(f"Skipping line due to ValueError: {line.strip()}")

with open('overlap_utr.txt', 'w') as overlap_outfile:
    overlap_outfile.writelines(overlap_lines)

with open('nonoverlap_utr.txt', 'w') as non_overlap_outfile:
    non_overlap_outfile.writelines(non_overlap_lines)

input_file = 'nonoverlap_utr.txt'
output_lines = []

with open(input_file, 'r') as infile:
    for line in infile:
        columns = line.strip().split('\t')
        if len(columns) >= 11:
            try:
                pA = int(columns[10])
                strand = columns[5]
                if strand == '+':
                    last = int(columns[2]) - 2000
                    first = int(columns[3])
                    cdsend = int(columns[1])
                    if pA < first:
                        output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{pA}\t{last}\t{cdsend}\n"
                    elif first < pA < last:
                        smaller = pA - first
                        larger = last - pA
                        if smaller <= larger:
                            output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{pA}\t{last}\t{cdsend}\n"
                        else:
                            output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{first}\t{pA}\t{cdsend}\n"
                    else:
                        output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{first}\t{pA}\t{cdsend}\n"

                elif strand == '-':
                    last = int(columns[1]) + 2000
                    first = int(columns[3])
                    cdsend = int(columns[2])
                    if pA > first:
                        output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{pA}\t{last}\t{cdsend}\n"
                    elif last < pA < first:
                        smaller = first - pA
                        larger = pA - last
                        if smaller <= larger:
                            output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{pA}\t{last}\t{cdsend}\n"
                        else:
                            output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{first}\t{pA}\t{cdsend}\n"
                    else:
                        output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{first}\t{pA}\t{cdsend}\n"

                else:
                    continue
                output_lines.append(output_line)

            except ValueError:
                print(f"Skipping line due to ValueError: {line.strip()}")

with open('processed_nonoverlap.tsv', 'w') as outfile:
    outfile.writelines(output_lines)

cmd = r"""awk '{print $5}' update_pas.txt | sort | uniq -c | awk -v id='1' '$1 != id {print $2}' > update_2+.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"""for i in $(cat update_2+.txt); do awk -v gene="$i" '$5 == gene {print $0}' update_pas.txt; done | sed 's/:/\t/g' > update_2+pas.txt"""
print('CMD:', cmd)
os.system(cmd)

data = pd.read_csv('update_2+pas.txt', delim_whitespace=True, header=None)
data.columns = range(data.shape[1])
gene_col = 4
value_col = 12
result = data.groupby(gene_col).apply(lambda x: x.nlargest(2, value_col)).reset_index(drop=True)
result.to_csv('update_2pas.txt', sep=' ', index=False, header=False)

cmd = r"""for i in $(cat update_2+.txt); do awk -v gene="$i" '$5 == gene {printf "%s ", $0}' update_2pas.txt; echo; done | sed 's/:/ /g' > update_2+pas.txt"""
print('CMD:', cmd)
os.system(cmd)

input_file = 'update_2+pas.txt'
overlap_lines = []
non_overlap_lines = []

with open(input_file, 'r') as infile:
    for line in infile:
        columns = line.strip().split()
        if len(columns) >= 25:
            try:
                pA1 = int(columns[10])
                pA2 = int(columns[24])
                strand = columns[5]

                if strand == '+':
                    last = int(columns[2]) - 2000
                    first = int(columns[3])
                elif strand == '-':
                    last = int(columns[1]) + 2000
                    first = int(columns[3])
                else:
                    continue
                if (abs(pA1 - first) <= 25 and abs(pA2 - last) <= 25) or \
                        (abs(pA1 - last) <= 25 and abs(pA2 - first) <= 25):
                    overlap_lines.append(line)
                else:
                    non_overlap_lines.append(line)

            except ValueError:
                print(f"Skipping line due to ValueError: {line.strip()}")

with open('overlap_2+utr.txt', 'w') as outfile:
    outfile.writelines(overlap_lines)

with open('nonoverlap_2+utr.txt', 'w') as outfile:
    outfile.writelines(non_overlap_lines)

input_file = 'nonoverlap_2+utr.txt'
output_lines = []

with open(input_file, 'r') as infile:
    for line in infile:
        columns = line.strip().split(' ')
        if len(columns) >= 25:
            try:
                pA1 = int(columns[10])
                pA2 = int(columns[24])
                strand = columns[5]

                smaller = min(pA1, pA2)
                larger = max(pA1, pA2)

                if strand == '+':
                    output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{smaller}\t{larger}\t{columns[1]}\n"
                elif strand == '-':
                    output_line = f"{columns[4]}\t{columns[0]}\t{strand}\t{larger}\t{smaller}\t{columns[2]}\n"
                else:
                    continue
                output_lines.append(output_line)

            except ValueError:
                print(f"Skipping line due to ValueError: {line.strip()}")

with open('processed_nonoverlap_2+.tsv', 'w') as outfile:
    outfile.writelines(output_lines)

cmd = r"""cat processed_nonoverlap.tsv processed_nonoverlap_2+.tsv | awk '{print $1}' > updateUTR_gene.txt"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"""for i in $(cat refUTR_gene.txt updateUTR_gene.txt | sort | uniq -u); do j=$(echo "\"$i\""); awk -v id="$j" '$1 == id {print $0}' """ + library_dir + """/apalyzer/refUTRraw_hg19.tsv; done > nonchangeUTR.tsv"""
print('CMD:', cmd)
os.system(cmd)

cmd = """head -1 """ + tools_dir + r"""/apalyzer/refUTRraw_hg19.tsv | sed 's/"//g' > updUTRraw_hg19.tsv"""
print('CMD:', cmd)
os.system(cmd)

cmd = r"""cat nonchangeUTR.tsv processed_newUTR.tsv processed_nonoverlap.tsv processed_nonoverlap_2+.tsv | sed 's/"//g' | sort -k2,2 -k4,4n >> updUTRraw_hg19.tsv"""
print('CMD:', cmd)
os.system(cmd)

numpy2ri.activate()
with localconverter(default_converter + numpy2ri.converter):
    robjects.r.source(library_dir+'/apalyzer_quant.R')
    result = \
        robjects.r.apalyzer_quant(out_dir)[0]
