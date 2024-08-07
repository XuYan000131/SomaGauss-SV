import os
from .config import minimap2,samtools,ms_dir,somatic_dir,hg38_fa,sniffles2,tr_bed


def bam_tag(script_tag,fq,new_ID,somatic_sample_dir,hg38_fa,stdout,stderr):
    os.system('/bin/bash {script} {a} {b} {c} {d} {minimap2} {samtools} 1 > {out} 2 > {err}'.format(
        script=script_tag, a=fq, b=new_ID, c=somatic_sample_dir, d=hg38_fa,minimap2=minimap2,samtools=samtools,out=stdout, err=stderr))


def Bulit_SH_Merge_and_call_SV(Sample):
    sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
    if not os.path.exists(sniffles2_sample_dir):
        os.makedirs(sniffles2_sample_dir)
    sv_vcf2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2.2.26_sniffles_v2.vcf' % Sample)
    bam_tag_dir_tumor = os.path.join(os.path.join(somatic_dir, 'tumor', Sample),
                                     '{DRR258593}_tumor_minimap2.2.26_sorted_tag.bam'.format(DRR258593=Sample))
    bam_tag_dir_blood = os.path.join(os.path.join(somatic_dir, 'blood', Sample),
                                     '{DRR258594}_blood_minimap2.2.26_sorted_tag.bam'.format(DRR258594=Sample))
    bam_merge_file = os.path.join(somatic_dir, 'Merge_bam', Sample)
    if not os.path.exists(bam_merge_file):
        os.makedirs(bam_merge_file)
    bam_merge_T_N = os.path.join(bam_merge_file, '%s_minimap2_sorted_T_N_merge.bam' % Sample)
    script_ms = os.path.join(ms_dir, Sample, '%s_bam_merge_sniffles.2nd.sh' % Sample)
    with open(script_ms, 'w') as out:
        out.write("#! /bin/bash" + '\n')
        out.write('''echo "$(date) 1. Start to merge bam: %s" ''' % Sample + '\n')
        out.write("{samtools} merge -@ 40 -h {bam_N} {out_bam} {bam_N} {bam_T}".format(
            samtools=samtools, bam_N=bam_tag_dir_blood, bam_T=bam_tag_dir_tumor, out_bam=bam_merge_T_N) + '\n')
        out.write('''/NAS/wg_fzt/software/samtools-1.9/samtools index -@ 25 %s \n''' % bam_merge_T_N)
        out.write('''echo "$(date) 1. Finish to merge bam: %s" ''' % Sample + '\n')
        out.write('''echo "$(date) 2. Start to sniffles2: %s" ''' % Sample + '\n')
        out.write(
            "{sniffles} -i {bam_sort} -v {vcf} --tandem-repeats {tr} -t 100 --minsupport 3 --mapq 50 --min-alignment-length 1000 "
            "--output-rnames --allow-overwrite --long-ins-length 100000 --cluster-binsize 500 --minsvlen 50 --reference {ref}".format(
                sniffles=sniffles2, bam_sort=bam_merge_T_N, vcf=sv_vcf2_T_N, tr=tr_bed, ref=hg38_fa) + '\n')
        out.write('''echo "$(date) 2. Finish to sniffles2: %s" ''' % Sample + '\n')



def ms_bam(script_ms_x,stdout,stderr):
    os.system( '/bin/bash {script} 1 > {out} 2 > {err} '.format(
        out=stdout, err=stderr, script=script_ms_x))


def change_vcf(path_vcf, work_dir,support_reads):
    f = path_vcf
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    vcfDf.columns = records[0]
    vcfDf['PRECISE'] = vcfDf['INFO'].apply(lambda x: x.split(";")[0])
    vcfDf['SVLen'] = vcfDf['INFO'].apply(
        lambda x: np.abs(int(x.split("SVLEN=")[-1].split(";")[0])) if re.search("SVLEN=", x) else 0)
    vcfDf['SVType'] = vcfDf['INFO'].apply(lambda x: x.split("SVTYPE=")[-1].split(";")[0])
    vcfDf['RNAMES'] = vcfDf['INFO'].apply(lambda x: x.split("RNAMES=")[-1].split(";")[0])
    ### Extract Tumor / Normal Support Reads
    vcfDf['TR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('tumor', r)], axis=1)
    vcfDf['NR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('blood', r)], axis=1)
    vcfDf['Support'] = vcfDf.apply(lambda x: x['TR'] + x['NR'], axis=1)
    vcfDf['TS'] = vcfDf['TR'].apply(lambda x: len(x))
    vcfDf['NS'] = vcfDf['NR'].apply(lambda x: len(x))
    vcfDf['TR'] = vcfDf.apply(lambda x: ",".join([r for r in x['RNAMES'].split(",") if re.search('tumor', r)]), axis=1)
    vcfDf['NR'] = vcfDf.apply(lambda x: ",".join([r for r in x['RNAMES'].split(",") if re.search('blood', r)]), axis=1)
    vcfDf['STDEV_POS'] = vcfDf['INFO'].apply(lambda x: x.split("STDEV_POS=")[-1].split(";")[0])
    # print(vcfDf)
    vcfDf['STDEV_LEN'] = vcfDf['INFO'].apply(
        lambda x: x.split("STDEV_LEN=")[-1].split(";")[0] if x.split("SVTYPE=")[-1].split(";")[0] not in ['BND'] else 0)
    vcfDf['STD'] = vcfDf.apply(lambda x: float(x['STDEV_POS']) + float(x['STDEV_LEN']), axis=1)
    vcfDf['POS2'] = vcfDf.apply(
        lambda x: int(x['POS']) + 1 if (x['SVType'] == 'INS') or (x['SVType'] == 'BND') else int(x['POS']) + int(
            x['SVLen']), axis=1)
    vcfDf2 = vcfDf
    vcfDf2 = vcfDf2[vcfDf2['NS'] == 0]
    vcfDf2 = vcfDf2[vcfDf2['TS'] >= support_reads] #support_reads = 3
    vcfDf4 = pd.DataFrame(vcfDf2,
                          columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    vcfDf5 = pd.DataFrame(vcfDf2,
                          columns=['#CHROM', 'POS', 'POS2', 'ID', 'SVType', 'SVLen', 'TS', 'NS', 'PRECISE', 'STD'])
    vcf_path = os.path.join(work_dir, os.path.basename(path_vcf).replace('.vcf', '.somaticSV.NS0_TS3.vcf'))
    vcfDf4['INFO'] = vcfDf4['INFO'].apply(lambda x: x.replace('SVTYPE=DUP', 'SVTYPE=INS'))
    vcfDf4.to_csv(vcf_path, index=False, header=False, sep='\t')
    bed_path = os.path.join(work_dir, os.path.basename(path_vcf).replace('.vcf', '.somaticSV.NS0_TS3.bed'))
    vcfDf5.to_csv(bed_path, index=False, header=False, sep='\t')
    ######################################  分成bed格式查看
    for SV_type in ['DEL','BND','INV','INS','DUP']:
        vcfDf5_SV_type = vcfDf5[vcfDf5['SVType'] == SV_type]
        bed_path_SV_type = os.path.join(work_dir, os.path.basename(path_vcf).replace('.vcf', '.somaticSV.NS0_TS3.{SV_type}.bed'.format(SV_type=SV_type)))
        vcfDf5_SV_type.to_csv(bed_path_SV_type, index=False, header=False, sep='\t')
    ######################################
    # vcfDf5_SV_type = vcfDf5[(vcfDf5['SVType'] == 'INS') | (vcfDf5['SVType'] == 'DUP')]
    # bed_path_SV_type = os.path.join(work_dir, os.path.basename(path_vcf).replace('.vcf','.somaticSV.NS0_TS3.{SV_type}.bed'.format(SV_type='INS')))
    # vcfDf5_SV_type.to_csv(bed_path_SV_type, index=False, header=False, sep='\t')
    ######################################
    os.system('head -n 70 {path_vcf} > {path_vcf}.head'.format(path_vcf=path_vcf))
    vcf_path_head = os.path.join(work_dir, os.path.basename(path_vcf).replace('.vcf', '.somaticSV.NS0_TS3.head.vcf'))
    os.system('cat {path_vcf}.head  {vcf_path} > {vcf_path_head}'.format(path_vcf=path_vcf,vcf_path=vcf_path, vcf_path_head=vcf_path_head))















