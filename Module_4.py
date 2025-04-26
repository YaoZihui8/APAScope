import glob, os, sys, time
import pandas as pd
import numpy as np
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects import StrVector
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects

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