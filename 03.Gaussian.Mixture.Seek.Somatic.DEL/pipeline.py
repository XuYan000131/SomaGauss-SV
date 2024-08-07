import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from .function import get_merge_bam_path,make_new_df_from_csv,get_tumor_support_reads_bam,get_blood_bam,CallSomaticP
from .config import ms_dir,POOL_core,bedtools,Sample





sample_name = Sample

## DEL does not need to construct CS sequences between the filtered abnormal to find

csv_path = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV.NS0_TS3.DEL.bed'.format(sampel_name=sample_name))
work_dir_000 = os.path.join(ms_dir, sample_name,'04.polish.DEL.GaussianMixture')

if not os.path.exists(work_dir_000):
    os.makedirs(work_dir_000)

tumor_bam_path, blood_bam_path = get_merge_bam_path(sample_name)
vcf_region_df = make_new_df_from_csv(csv_path, sample_name) ## Specific pairs of DELs constructed for the situation

ParamList = [(work_dir_000,tumor_bam_path,blood_bam_path,vcf_region_df, index) for index in vcf_region_df.index]
print('Start extracting the bam file for tumor')
P = Pool(POOL_core)
P.map(get_tumor_support_reads_bam, ParamList) #Extracting the bam and fq files for tumor support sv returns no value.
P.close()
del P
print('Finish extracting the bam file for tumor')
print('Start extracting bloods bam file')
P = Pool(POOL_core)
P.map(get_blood_bam, ParamList) #提取blood的bam文件 没有返回值
P.close()
del P
print('Finished extracting bloods bam file')
print('Start bringing in a Gaussian distribution model to check if there is a true')
P = Pool(POOL_core)
RecordList = P.map(CallSomaticP, ParamList)
P.close()
del P
print('output result')
JudgementDf = pd.DataFrame.from_records(RecordList)
JudgementDf.columns = ['dir_name', 'JUDGEMENT']
JudgementDf.index = np.array(JudgementDf['dir_name'])
ComDf = pd.concat([vcf_region_df, JudgementDf], axis=1) #Merge two dataframes with dir_name as index
print('Finished extracting bloods bam file')
ComDf.to_csv(os.path.join(work_dir_000, '%s_GaussianMixture_Judgement.csv' % (sample_name)), sep='\t', index=False)
ComDf = ComDf.sort_values(by=['chr', 'start'], ascending=True)
ComDf_Ture = ComDf[ComDf['JUDGEMENT'] == 'Ture_somatic_DEL']
ComDf_bed = ComDf_Ture[['chr', 'start', 'end', 'sv_id']]
ComDf_bed.to_csv(os.path.join(work_dir_000, '%s_GaussianMixture_Judgement.bed' % (sample_name)), sep='\t', index=False,header=False)  ### 只提取判断为真的SV 情况
os.system('''{bedtools} cluster -d 200 -i {bed}   > {bed_cluster}'''.format(bedtools=bedtools,bed =os.path.join(work_dir_000, '%s_GaussianMixture_Judgement.bed' % (sample_name)),
                                                                            bed_cluster=os.path.join(work_dir_000, '%s_GaussianMixture_Judgement_cluster.bed' % (sample_name))))





