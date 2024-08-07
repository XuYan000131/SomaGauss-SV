import os


samtools = "samtools" # Version 1.9
bedtools = "bedtools" # Version 2.30.0
sniffles2 = "sniffles" #Version 2.0.6
minimap2 = 'minimap2' #Version 2.26
seqkit = 'seqkit' #Version 2.4.0
tr_bed = "human_GRCh38_no_alt_analysis_set.trf.bed" #from sniffels
hg38_fa = "hg38_mainChr.fa"
Repeat_masked_bed = "mainChr_hg38_UCSC.RepeatMasked.Unselected.bed"
hg38_fa_bed = "hg38_mainChr.fa.bed"

straglr = 'straglr.py' #Version 1.4.1
conda_activate = '/home/bin/activate'
straglr_env = '/home/miniconda3/envs/straglr'




work_dir = "work_dir"
POOL_core = 4
support_reads = 3
Sample = 'HCC1395'
tumor_fq = 'tumor_fq_path'
blood_fq = 'blood_fq_path'

program_dir = os.path.join(work_dir, '01.program')
ms_dir = os.path.join(work_dir, '11.minimap2_sniffles')
coverage_dir = os.path.join(work_dir,'03.coverage')
somatic_dir = os.path.join(work_dir, "04.Somatic_minimap2_2.4")
work_dir_count = os.path.join(work_dir, "05.base_count")
