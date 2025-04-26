import glob, os, sys, time
import pandas as pd
import numpy as np
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects import StrVector
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects

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


