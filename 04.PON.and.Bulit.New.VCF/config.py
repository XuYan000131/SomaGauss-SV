
import os
import re
import pandas as pd
from multiprocessing import Pool


work_dir = "work_dir"
ms_dir = os.path.join(work_dir, '11.minimap2_sniffles_DeBreak')
Sample = 'HCC1395'


jasmine = "jasmine" #Version 1.1.5
hg38_fa = "hg38_mainChr.fa"

XieZhi_405_PON_vcf = 'Chinese_405_SV_geno.seq.new.20240712.vcf'
