import glob, os, sys, time
import pandas as pd
import numpy as np
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects import StrVector
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects

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