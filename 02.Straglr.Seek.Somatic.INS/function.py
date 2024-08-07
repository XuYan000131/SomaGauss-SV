
import os
from .config import ms_dir,somatic_dir
from .config import hg38_fa,Repeat_masked_bed,straglr,straglr_env,conda_activate,bedtools
import pandas as pd


def Cat_RepeatMasked_IN_OR_NOT_IN_funcatin(sample_name):
    # path_vcf = '/NAS/wg_liuxy/01.Pancaner/Iris_test/{OSCC5}/{OSCC5}_iris_out.vcf'.format(OSCC5=sample_name)
    ## DEL 不用构建CS序列之间在筛选abnormal的里面去找
    csv_path_INS = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV.NS0_TS3.INS.bed'.format(sampel_name=sample_name))
    csv_path_DUP = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV.NS0_TS3.DUP.bed'.format(sampel_name=sample_name))
    csv_path_DUP_df = pd.read_csv(csv_path_DUP,sep='\t',header=None)
    csv_path_DUP_df[2] = csv_path_DUP_df[1] + 1
    csv_path_DUP_2 = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_DUP_1_bp.bed'.format(sampel_name=sample_name))
    csv_path_DUP_df.to_csv(csv_path_DUP_2,sep='\t',header=None,index=None)
    csv_path = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP.bed'.format(sampel_name=sample_name))
    vcf_file = os.path.join(ms_dir,sample_name,'sniffles2','{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV.NS0_TS3.head.vcf'.format(sampel_name=sample_name))
    cmd_cat_INS_and_DUP = "cat {csv_path_INS} {csv_path_DUP} > {csv_path}".format(csv_path_INS=csv_path_INS, csv_path_DUP=csv_path_DUP_2, csv_path=csv_path)
    os.system(cmd_cat_INS_and_DUP)
    work_dir_INS = os.path.join(ms_dir, sample_name,'05.INS.SIX.KINDS')
    if not os.path.exists(work_dir_INS):
        os.makedirs(work_dir_INS)
    Intersect_bed_path = os.path.join(work_dir_INS,
                                      '{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP_with_UCSC_RepeatMasked.bed'.format(
                                          sampel_name=sample_name))
    NO_Intersect_bed_path = os.path.join(work_dir_INS,
                                         '{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP_NONE_UCSC_RepeatMasked.bed'.format(
                                             sampel_name=sample_name))
    cmd_bedtools_inter = "{bedtools} intersect -a {csv_path} -b {Repeat_masked_bed} -wa -wb > {Intersect_bed_path}".format(
        bedtools=bedtools, csv_path=csv_path, Repeat_masked_bed=Repeat_masked_bed,
        Intersect_bed_path=Intersect_bed_path)
    cmd_bedtools_NO_inter = "{bedtools} intersect -a {csv_path} -b {Repeat_masked_bed} -wa -wb -v > {NO_Intersect_bed_path}".format(
        bedtools=bedtools, csv_path=csv_path, Repeat_masked_bed=Repeat_masked_bed,
        NO_Intersect_bed_path=NO_Intersect_bed_path)
    os.system(cmd_bedtools_inter)  # 判断是否和repeat有交集
    os.system(cmd_bedtools_NO_inter)  # 判断是否和repeat有交集
    ALL_Intersect_bed_path = os.path.join(work_dir_INS,'{sampel_name}_T_N_merge_minimap2.2.26_sniffles_v2.somaticSV_INS_and_DUP_all_UCSC_RepeatMasked.bed'.format(sampel_name=sample_name))
    os.system("cat {Intersect_bed_path} {NO_Intersect_bed_path} > {ALL_Intersect_bed_path}".format(Intersect_bed_path=Intersect_bed_path, NO_Intersect_bed_path=NO_Intersect_bed_path, ALL_Intersect_bed_path=ALL_Intersect_bed_path))



def make_new_df_from_Intersect_bed(csv_path, sample_name):
    csv_df = pd.read_csv(csv_path, sep='\t')
    csv_df.columns = ['chr', 'start', 'end', 'svid','svtype', 'svlen', 'support_number','number', 'PRECISE', 'STD',
                      'chr_repeat', 'start_repeat', 'end_repeat','strand','sub_repeat_type','repeat_type','number_1',
                      'number_2','number_3']
    csv_df = csv_df.fillna('NA')
    csv_df_DEL = csv_df.loc[(csv_df['svtype'] == 'INS')|(csv_df['svtype'] == 'DUP')]
    csv_df_DEL_region = pd.DataFrame()
    csv_df_DEL_region['chr'] = csv_df_DEL['chr']
    csv_df_DEL_region['start'] = csv_df_DEL.apply(lambda x: x['start']-500 if x['start_repeat']=='NA' else x['start_repeat']-500, axis=1)
    csv_df_DEL_region['end'] = csv_df_DEL.apply(lambda x: x['end']+500 if x['end_repeat']=='NA'  else x['end_repeat']+500, axis=1)
    csv_df_DEL_region['start'] = csv_df_DEL_region['start'].astype(int)
    csv_df_DEL_region['end'] = csv_df_DEL_region['end'].astype(int)
    csv_df_DEL_region['sv_id'] = csv_df_DEL['svid']
    csv_df_DEL_region['sample_name'] = sample_name
    csv_df_DEL_region = csv_df_DEL_region.sort_values(by=['chr', 'start', 'end'])
    return csv_df_DEL_region



def get_merge_bam_path_tumor_and_blood(sample_name):
    tumor_bam = os.path.join(os.path.join(somatic_dir, 'tumor', sample_name),
                                     '{DRR258593}_tumor_minimap2.2.26_sorted_tag.bam'.format(DRR258593=sample_name))
    blood_bam = os.path.join(os.path.join(somatic_dir, 'blood', sample_name),
                                     '{DRR258594}_blood_minimap2.2.26_sorted_tag.bam'.format(DRR258594=sample_name))
    return (tumor_bam,blood_bam)



def straglr_function(work_dir,sample_name,vcf_region_bed_2_path,blood_bam_path,tumor_bam_path):
    straglr_out_dir_blood = os.path.join(work_dir, '01.blood.straglr_out/{sampel}_blood.straglr'.format(sampel=sample_name))
    tmp_dir_blood = os.path.join(work_dir, '01.blood.straglr_out', 'tmp')
    if not os.path.exists(tmp_dir_blood):
        os.makedirs(tmp_dir_blood)
    cmd_straglr_blood = 'source {conda_activate} {straglr_env} && {straglr} {bam} {hg38_fa} {straglr_out_dir} ' \
                        '--regions {vcf_region_bed_2_path} --genotype_in_size --min_support 1 --max_str_len 50 --nprocs 200 --min_ins_size 50 --max_num_clusters {cluster}  --tmpdir {tmp_dir}'.format(straglr=straglr,conda_activate=conda_activate,
    straglr_env=straglr_env,bam=blood_bam_path, hg38_fa=hg38_fa, straglr_out_dir=straglr_out_dir_blood, vcf_region_bed_2_path=vcf_region_bed_2_path, cluster=9, tmp_dir=tmp_dir_blood)
    os.system(cmd_straglr_blood)
    straglr_out_dir_tumor = os.path.join(work_dir, '02.tumor.straglr_out/{sampel}_tumor.straglr'.format(sampel=sample_name))
    tmp_dir_tumor = os.path.join(work_dir, '02.tumor.straglr_out', 'tmp')
    if not os.path.exists(tmp_dir_tumor):
        os.makedirs(tmp_dir_tumor)
    cmd_straglr_tumor = 'source {conda_activate} {straglr_env} && {straglr} {bam} {hg38_fa} {straglr_out_dir} ' \
                        '--regions {vcf_region_bed_2_path} --genotype_in_size --min_support 1 --max_str_len 100 --nprocs 200 --min_ins_size 50 --max_num_clusters {cluster}  --tmpdir {tmp_dir}'.format(straglr=straglr,conda_activate=conda_activate,
                        straglr_env=straglr_env,bam=tumor_bam_path, hg38_fa=hg38_fa, straglr_out_dir=straglr_out_dir_tumor, vcf_region_bed_2_path=vcf_region_bed_2_path, cluster=9, tmp_dir=tmp_dir_tumor)
    os.system(cmd_straglr_tumor)




def Filter_somatic(straglr_out):
    straglr_out['counts'] = straglr_out['genotype'].apply(lambda x:len(x.split(';')))
    straglr_out = straglr_out[straglr_out['counts']>1]
    straglr_df = straglr_out.groupby(['#chrom','start','end','locus'])['read_name'].agg(','.join).reset_index()
    straglr_df['TR'] = straglr_df.apply(lambda x: [r for r in x['read_name'].split(",") if re.search('tumor',r)], axis=1)
    straglr_df['NR'] = straglr_df.apply(lambda x: [r for r in x['read_name'].split(",") if re.search('blood',r)], axis=1)
    straglr_df['TS'] = straglr_df['TR'].apply(lambda x: len(x))
    straglr_df['NS'] = straglr_df['NR'].apply(lambda x: len(x))
    straglr_df['TR'] = straglr_df.apply(lambda x: ",".join([r for r in x['TR']]), axis=1)
    straglr_df['NR'] = straglr_df.apply(lambda x: ",".join([r for r in x['NR']]), axis=1)
    ### Filter somatic
    straglr_df['somatic'] = straglr_df.apply(lambda x: 'somatic' if x['TS']>=3 and x['NS']==0 else 'germline', axis=1)
    return(straglr_df)



