


import os
import re
import pandas as pd
from multiprocessing import Pool







def change_vcf(path_vcf):
    f = path_vcf
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    vcfDf.columns = records[0]
    return vcfDf




#run sh
def ms_bam(script_ms_x,stdout,stderr):
    os.system( '/bin/bash {script} 1 > {out} 2 > {err} '.format(
        out=stdout, err=stderr, script=script_ms_x))

