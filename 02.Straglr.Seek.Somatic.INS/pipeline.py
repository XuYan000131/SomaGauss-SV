

import os
import pandas as pd
from .function import Cat_RepeatMasked_IN_OR_NOT_IN_funcatin,get_merge_bam_path_tumor_and_blood,make_new_df_from_Intersect_bed,Filter_somatic,straglr_function
from .config import ms_dir,bedtools,Sample



#### Starting with a little more than one step, intersecting with Repeatmasker, and then catting it together ##
Cat_RepeatMasked_IN_OR_NOT_IN_funcatin(Sample)

#####  Straglr tumorbam and bloodbam separately and then filter somatic INS
sample_name = Sample
work_dir_straglr = os.path.join(ms_dir, sample_name,'05.INS.SIX.KINDS')
all_VCF_region_df= pd.read_csv(os.path.join(work_dir_straglr, 'Finally_{sample_name}_all_VCF_region_df.bed'.format(sample_name=sample_name)),header = 0,sep='\t')
ALL_Intersect_bed_path = os.path.join(work_dir_straglr,'{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP_all_UCSC_RepeatMasked.bed'.format(sampel_name=sample_name))
ALL_Intersect_bed_path_df = pd.read_csv(ALL_Intersect_bed_path,header = None,sep='\t')
ALL_Intersect_bed_path_bed_path = os.path.join(work_dir_straglr,'{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP_all_UCSC_RepeatMasked.chr.star.end.id.bed'.format(sampel_name=sample_name))
ALL_Intersect_bed_path_df[[0,1,2,3]].to_csv(ALL_Intersect_bed_path_bed_path, index=False,sep='\t',header=False)
tumor_bam_path, blood_bam_path = get_merge_bam_path_tumor_and_blood(sample_name) #, blood_bam_path
vcf_region_bed = make_new_df_from_Intersect_bed(ALL_Intersect_bed_path, sample_name) ## Specific to the case of INS constructs
vcf_region_bed_2 = vcf_region_bed[vcf_region_bed['sv_id'].isin(all_VCF_region_df['sv_id'])]
vcf_region_bed_2_path = os.path.join(work_dir_straglr,'{sampel_name}_INS_and_DUP_all_UCSC_RepeatMasked.bed'.format(sampel_name=sample_name))
vcf_region_bed_2.to_csv(vcf_region_bed_2_path, index=False,sep='\t',header=False)
straglr_function(work_dir_straglr,sample_name,vcf_region_bed_2_path,blood_bam_path,tumor_bam_path) ###Distribution of blood and tumor for straglr analysis
tumor_straglr_out_tsv = os.path.join(work_dir_straglr, '02.tumor.straglr_out/{sampel}_tumor.straglr.tsv'.format(sampel=sample_name))
blood_straglr_out_tsv = os.path.join(work_dir_straglr, '01.blood.straglr_out/{sampel}_blood.straglr.tsv'.format(sampel=sample_name))
tumor_straglr_out_tsv_df = pd.read_csv(tumor_straglr_out_tsv, header=1, sep='\t')
blood_straglr_out_tsv_df = pd.read_csv(blood_straglr_out_tsv, header=1, sep='\t')
all_straglr_out_tsv_df = pd.concat([tumor_straglr_out_tsv_df, blood_straglr_out_tsv_df], axis=0)
all_straglr_out_tsv_df = all_straglr_out_tsv_df[all_straglr_out_tsv_df['read_status'] != 'skipped (not_spanning)']#read_status
all_straglr_out_tsv_df = all_straglr_out_tsv_df.sort_values(by=['locus','read_name'])
straglr_df = Filter_somatic(all_straglr_out_tsv_df)
straglr_df.sort_values(by=['#chrom','start']).to_csv(os.path.join(work_dir_straglr, '03.straglr_out_{sampel}_df.tsv'.format(sampel=sample_name)), index=False, sep='\t')
straglr_df[['#chrom','start','end','somatic']].sort_values(by=['#chrom','start']).to_csv(os.path.join(work_dir_straglr, '04.straglr_out_{sampel}_df.bed'.format(sampel=sample_name)), index=False, sep='\t',header=False)
# all_VCF_region_df.to_csv(os.path.join(work_dir, '05.Finally_{sample_name}_all_VCF_region_df.bed'.format(sample_name=sample_name)), index=False, sep='\t',header=False)
somatic_straglr_path = os.path.join(work_dir_straglr, '06.straglr_somatic_out_{sampel}_intersect.bed'.format(sampel=sample_name))
cmd_intersect = '{bedtools} intersect -a {ALL_Intersect_bed_path_bed_path} -b {straglr_bed} -wa -wb > {somatic_straglr_path}'\
    .format(bedtools=bedtools,ALL_Intersect_bed_path_bed_path=ALL_Intersect_bed_path_bed_path,straglr_bed=os.path.join(work_dir_straglr, '04.straglr_out_{sampel}_df.bed'.format(sampel=sample_name)),
            somatic_straglr_path=somatic_straglr_path)
os.system(cmd_intersect)
somatic_straglr_df = pd.read_csv(somatic_straglr_path, header=None, sep='\t')
somatic_straglr_df_3_7 = somatic_straglr_df[[3,7]]
somatic_straglr_df_3_7.columns = ['sv_id', 'somatic']
all_VCF_region_df_somatic_judge = pd.merge(all_VCF_region_df, somatic_straglr_df_3_7, how='left', left_on='sv_id', right_on='sv_id')
all_VCF_region_df_somatic_judge.to_csv(os.path.join(work_dir_straglr, '07.Finally_{sample_name}_all_VCF_region_df_somatic_judge.bed'.format(sample_name=sample_name)), index=False, sep='\t')
all_VCF_region_df_somatic_judge.to_csv(os.path.join(ms_dir, '05.INS.SIX.KINDS','08.Finally_{sample_name}_all_VCF_region_df_somatic_judge.tsv'.format(sample_name=sample_name)), index=False, sep='\t')

###################################################################################################################################################################################################################################################################################################
#############################                                                                                                                                                                          #######################################################################
#####################################################################################################################################################################################################################################################################################################
###  Basic statistics

all_sample_classify = pd.DataFrame()

work_dir = os.path.join(ms_dir, sample_name,'05.INS.SIX.KINDS')
all_VCF_region_df_somatic_judge = pd.read_csv(os.path.join(ms_dir, '05.INS.SIX.KINDS','08.Finally_{sample_name}_all_VCF_region_df_somatic_judge.tsv'.format(sample_name=sample_name)), sep='\t',header=0)
all_VCF_region_df_somatic_judge.fillna('Unknow', inplace=True)
all_VCF_region_df_somatic_judge_1 = all_VCF_region_df_somatic_judge[all_VCF_region_df_somatic_judge['somatic'] != 'germline']
all_VCF_region_df_somatic_judge_1.to_csv(os.path.join(work_dir, '08.Finally_{sample_name}_all_VCF_region_df_somatic_judge_1.tsv'.format(sample_name=sample_name)), index=False, sep='\t')
all_VCF_region_df_somatic_judge_1[sample_name] = all_VCF_region_df_somatic_judge_1['classify']
all_sample_classify = all_sample_classify.append(all_VCF_region_df_somatic_judge_1[sample_name].value_counts())




all_sample_classify.to_csv(os.path.join(ms_dir, '05.INS.SIX.KINDS','somatic_all_sample_INS_DUP_classify_count.csv'), index=True)


# Stacked Histogram Statistics
Y_list = all_sample_classify.columns
all_sample_classify['sample'] = all_sample_classify.index
all_sample_classify.plot.bar(stacked=True,x = 'sample',y = Y_list,figsize=(10, 6))
png_file = os.path.join(ms_dir, '05.INS.SIX.KINDS','somatic_all_sample_INS_DUP_classify_count.png')
plt.savefig(png_file, dpi=300, bbox_inches='tight')

