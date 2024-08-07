

import os
import re
import pandas as pd
from multiprocessing import Pool
from .function import change_vcf,ms_bam
from .config import work_dir,ms_dir,jasmine,hg38_fa,XieZhi_405_PON_vcf,Sample



# Build a new VCF (vcf build after filtering)

INS_tsv_path = os.path.join(ms_dir, Sample,'05.INS.SIX.KINDS','08.Finally_{sample_name}_all_VCF_region_df_somatic_judge_1.tsv'.format(sample_name=Sample))
INS_df = pd.read_csv(INS_tsv_path,sep='\t',header=0)
DEL_bed_path = os.path.join(ms_dir, Sample,'04.polish.DEL.GaussianMixture', '%s_GaussianMixture_Judgement.bed' % (Sample))
DEL_df = pd.read_csv(DEL_bed_path, sep='\t', header=None,names=['chr','start','sv_id'])
sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
head_vcf_path = os.path.join(sniffles2_sample_dir,'%s_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV.NS0_TS3_STD25.head.vcf'% (Sample))
path_vcf= os.path.join(sniffles2_sample_dir,'%s_T_N_merge_minimap2.2.26_sniffles_v2.vcf'% (Sample))
Filter_DEL_INS_list =DEL_df['sv_id'].to_list()+INS_df['sv_id'].to_list()
vcfDf = change_vcf(head_vcf_path)
vcfDf_filter = vcfDf[(vcfDf['ID'].isin(Filter_DEL_INS_list))|(vcfDf['ID'].str.contains('BND'))|(vcfDf['ID'].str.contains('INV'))|(vcfDf['ID'].str.contains('TRA'))]
vcfDf_filter.to_csv(os.path.join(sniffles2_sample_dir,'somaticSV.%s_T_N_merge_minimap2.2.26_sniffles_v2.Finally.bed'% (Sample)),sep='\t',header=None,index=None)
vcf_path_head = os.path.join(sniffles2_sample_dir,'somaticSV.%s_T_N_merge_minimap2.2.26_sniffles_v2.Finally.vcf'% (Sample))
os.system('cat {path_vcf}.head  {vcf_path} > {vcf_path_head}'.format(path_vcf=path_vcf, vcf_path=os.path.join(sniffles2_sample_dir,'somaticSV.%s_T_N_merge_minimap2.2.26_sniffles_v2.Finally.bed'% (Sample)),
                                                                     vcf_path_head=vcf_path_head))







sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
vcf_path_head = os.path.join(sniffles2_sample_dir,
                             'somaticSV.%s_T_N_merge_minimap2.2.26_sniffles_v2.Finally.vcf' % (Sample))
Jasmine_dir = os.path.join(sniffles2_sample_dir,'02.Jasmine_dir')
if not os.path.exists(Jasmine_dir):
    os.makedirs(Jasmine_dir)
script_ms = os.path.join(Jasmine_dir,'%s.Jasmine.sh'% (Sample))
file_list  = os.path.join(Jasmine_dir,'%s.Jasmine.file.list.txt'% (Sample))
out_file =  os.path.join(Jasmine_dir, '%s.Jasmine.reslut.vcf'% (Sample))
with open(file_list, 'w') as out:
    out.write(vcf_path_head + "\n")
    out.write(XieZhi_405_PON_vcf + "\n")
with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write('''echo "$(date) 1. Start to jasmine : %s" ''' % Sample + '\n')
    cmd_jasmine = 'source  /home/bin/activate /home/miniconda3/envs/jasmine_env && {jasmine} threads=30 min_seq_id=0.9 k_jaccard=9 --dup_to_ins ignore_strand genome_file={ref} max_dist=500 file_list={sample_vcf} out_file={out_vcf}'.format(
        jasmine=jasmine, sample_vcf=file_list, out_vcf=out_file,ref=hg38_fa)
    out.write(cmd_jasmine + '\n')
    out.write('''echo "$(date) 1. Finish to jasmine : %s" ''' % Sample + '\n')




sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
Jasmine_dir = os.path.join(sniffles2_sample_dir, '02.Jasmine_dir')
script_ms = os.path.join(Jasmine_dir, '%s.Jasmine.sh' % (Sample))
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
ms_bam(script_ms, stdout, stderr)





Jasmine_count = pd.DataFrame()

sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
Jasmine_dir = os.path.join(sniffles2_sample_dir, '02.Jasmine_dir')
reslut_vcf = os.path.join(Jasmine_dir, '%s.Jasmine.reslut.vcf' % (Sample))
reslut_vcf_PON = os.path.join(Jasmine_dir, '%s.PON.reslut.vcf' % (Sample))
Cmd = "grep -E '#|Sniffles2' {reslut_vcf} > {reslut_vcf_PON}".format(reslut_vcf=reslut_vcf,reslut_vcf_PON=reslut_vcf_PON)
# os.system(Cmd)
f = reslut_vcf_PON
with open(f, 'r') as fin:
    records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
vcfDf = pd.DataFrame.from_records(records[1:])
columns_list = records[0]
vcfDf.columns = records[0]
vcfDf['SUPP_number'] = vcfDf.apply(lambda x:x['INFO'].split('SUPP=')[1].split(';')[0],axis =1)
vcfDf['ID'] = vcfDf.apply(lambda x:x['ID'].split('_')[1],axis =1)
Jasmine_count[Sample] = vcfDf['SUPP_number'].value_counts()
vcfDf_After_pon = vcfDf[vcfDf['SUPP_number'] == '1']
reslut_bed_PON = os.path.join(Jasmine_dir, 'somatic.%s.reslut.AF.PON.20240712.bed' % (Sample))
reslut_vcf_PON_head_txt = reslut_vcf_PON.replace('.vcf','.txt')
CMD_head = 'head -n 84 {reslut_vcf_PON} > {reslut_vcf_PON_head_txt}'.format(reslut_vcf_PON=reslut_vcf_PON,reslut_vcf_PON_head_txt=reslut_vcf_PON_head_txt)
os.system(CMD_head)
vcfDf_After_pon[columns_list].to_csv(reslut_bed_PON,sep='\t',header=None,index=False)
reslut_vcf_PON_2 = reslut_bed_PON.replace('.bed','.vcf')
CMD_cat = 'cat {reslut_vcf_PON_head_txt} {reslut_bed_PON} > {reslut_vcf_PON_2}'.format(reslut_bed_PON=reslut_bed_PON,reslut_vcf_PON_2=reslut_vcf_PON_2,
                                                                            reslut_vcf_PON_head_txt=reslut_vcf_PON_head_txt)
os.system(CMD_cat)


Jasmine_count.T.to_csv(os.path.join(ms_dir,'AF.PON.SV.number.count.20240712.txt'),sep='\t',header=True,index=True)

