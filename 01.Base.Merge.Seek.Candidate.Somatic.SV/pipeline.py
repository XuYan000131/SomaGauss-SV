
import os
from multiprocessing import Pool
from .function import bam_tag,Bulit_SH_Merge_and_call_SV,ms_bam,change_vcf
from .config import program_dir,ms_dir,coverage_dir,somatic_dir,POOL_core,hg38_fa
from .config import Sample,tumor_fq,blood_fq,support_reads

dir_list = [program_dir, ms_dir, somatic_dir,coverage_dir]
for dir in dir_list:
    if not os.path.exists(dir):
        os.mkdir(dir)


## Basic filtering and comparison of fastq
pools = Pool(POOL_core)

tag_list = ['tumor','blood']

for tag in tag_list:
    # sampletype_dir = os.path.join(data_dir, sampletype)
    if tag == 'tumor':
        fq = tumor_fq
    if tag == 'blood':
        fq = blood_fq
    script_tag = os.path.join(program_dir, "ONT_clean_fq_minmap2.2.26_%s.sh" % tag)  ## 在 01.Base.Merge.Seek.Candidate.Somatic.SV
    somatic_sample_dir = os.path.join(somatic_dir, tag, Sample)
    if not os.path.exists(somatic_sample_dir):
        os.makedirs(somatic_sample_dir)
    new_ID = '%s_%s' % (Sample, tag)
    stdout = os.path.join(somatic_sample_dir, '%s_add_tag.o' % new_ID)
    stderr = os.path.join(somatic_sample_dir, '%s_add_tag.e' % new_ID)
    pools.apply_async(bam_tag, args=(script_tag,fq,new_ID,somatic_sample_dir,hg38_fa,stdout,stderr))


pools.close()
pools.join()
del pools


## Building SH Scripts
Bulit_SH_Merge_and_call_SV(Sample)



## Run  sh script
script_ms = os.path.join(ms_dir,Sample, '%s_bam_merge_sniffles.2nd.sh' % Sample)
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms,stdout,stderr)



#############################################################################################################################################################################
############################    Use the script to filter the results of sniffles to get candidate somatic SVs.
############################################################################################################################################################################


sniffles2_sample_dir = os.path.join(ms_dir, Sample, 'sniffles2')
sv_vcf2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2.2.26_sniffles_v2.vcf' % Sample)
# pools.apply_async(change_vcf, args=(sv_vcf2_T_N, sniffles2_sample_dir,3)) #support_reads =3
change_vcf(sv_vcf2_T_N, sniffles2_sample_dir, support_reads)
















