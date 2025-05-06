import os
from sklearn.mixture import GaussianMixture
import scipy.stats as stats
import pandas as pd
import numpy as np
import pysam
import random
import matplotlib.pyplot as plt
from .config import samtools,somatic_dir

#############################################################################################################################################################################
############################  Filtering DEL using Gaussian distribution models
############################################################################################################################################################################


def PNG(Newdf,New_x,New_labels,New_pie,means,stds,length_list,labels,png_path):
    New_readsLabel = Newdf['T_or_N'].values
    normal_reads = New_x[np.where(New_readsLabel == 'N')[0]]
    normal_reads = normal_reads.reshape((normal_reads.shape[0],))
    likelihoods = []
    for mu, std in zip(means[np.unique(New_labels)], stds[np.unique(New_labels)]):
        likelihoods.append(np.exp(-0.5 * ((normal_reads - mu) / std) ** 2) + 1e-10)
    likelihoods = np.vstack(likelihoods).T
    # Calculate the posterior probability
    posteriors = (likelihoods * New_pie).sum(axis=0) / (likelihoods * New_pie).sum()
    x = np.linspace(min(length_list) - 20, max(length_list) + 20, 100)
    labels_list = set(labels.tolist())
    # Iterate through the list of means and variances and plot each normal distribution curve
    n = 0
    #Setting the frame size
    plt.figure(figsize=(10, 5))
    for mean, std in zip(means[np.unique(New_labels)], stds[np.unique(New_labels)]):
        if std <= 0.001:
            print('std is too small')
        else:
            y = (1 / (np.sqrt(2 * np.pi) * std)) * np.exp(-0.5 * ((x - mean) / std) ** 2)
            plt.plot(x, y, label=f"label= {n}")
            plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        n = n + 1
    plt.scatter(Newdf['sv_len'], np.ones(Newdf.shape[0]) * (-0.25),
                c=Newdf['T_or_N'].apply(lambda x: 'b' if x == 'N' else 'r'), alpha=0.5, s=10)
    plt.legend()
    plt.xlabel('X')
    plt.ylabel('Probability Density')
    plt.savefig(png_path)
    plt.close()
    # return(np.vstack([means[np.unique(New_labels)], stds[np.unique(New_labels)], AimCluster, posteriors]).T)



def CallSomaticP(ParamList):
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    print(index)
    return_relsult = 'return'
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    region = line_x['region']
    blood_bed_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_blood.csv'.format(sampel_name_x=dir_name))
    tumor_bed_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_tumor.csv'.format(sampel_name_x=dir_name))
    GaussianMixture_dir = os.path.join(wok_dir, dir_name, 'GaussianMixture')
    if not os.path.exists(GaussianMixture_dir):
        os.makedirs(GaussianMixture_dir)
    df_blood = pd.read_csv(blood_bed_out, sep='\t', header=0)
    df_tumor = pd.read_csv(tumor_bed_out, sep='\t', header=0)
    name_smple = dir_name
    df_all_path = os.path.join(GaussianMixture_dir, '%s_all.csv' % (name_smple))
    png_path = os.path.join(GaussianMixture_dir, 'GMM.%s.png' % (name_smple))
    if df_tumor.shape[0] < 3:
        return_relsult = 'False_somatic_DEL.tumor_support_reads_less_3.1'
        # print('False_somatic_DEL.tumor_support_reads_less_3.1')
    elif df_blood.shape[0] == 0:
        if df_tumor.shape[0] >= 3:
            return_relsult = 'Ture_somatic_DEL'
            # print('Ture_somatic_DEL')
        else:
            return_relsult = 'False_somatic_DEL.tumor_support_reads_less_3.0'
            # print('False_somatic_DEL.tumor_support_reads_less_3.0')
            # return ((dir_name, 'False_somatic_DEL.tumor_support_reads_less_3.0'))
    if return_relsult == 'return' : #or return_relsult == 'Ture_somatic_DEL_1'
        df_all = pd.concat([df_blood, df_tumor], axis=0)
        df_all['sv_len'] = df_all['sv_len'].apply(lambda x: x + random.uniform(0.001, 0.0001))
        length_list = df_all['sv_len'].values
        CountPool = length_list
        x_1111 = CountPool.reshape((len(CountPool), 1))
        # print('x_1111.shape', x_1111.shape)
        best_gmm, labels = SelectK(x_1111) ## Determining the optimal grouping
        # print('best_gmm.n_components', best_gmm.n_components)
        df_all['label'] = labels
        df_all['T_or_N'] = df_all.apply(lambda x: 'N' if 'blood' in x['query_name'] else 'T', axis=1)
        df_all.to_csv(df_all_path, sep='\t', index=False)
        means = best_gmm.means_.reshape((best_gmm.n_components,))
        covs = best_gmm.covariances_
        stds = np.array([np.sqrt(covs[i]) for i in range(best_gmm.n_components)]).reshape((best_gmm.n_components,))
        # df_all['label'] = labels
        cluster, count = np.unique(labels, return_counts=True)
        outlinerT = 3
        AimCluster = cluster[np.where(count >= outlinerT)[0]] ##  Filter groups with reads greater than three
        New_x = x_1111[np.where(np.in1d(labels, AimCluster))[0]]
        Newdf = df_all[df_all['label'].isin(AimCluster)]
        if Newdf.shape[0] == 0:
            return_relsult = 'False_somatic_DEL.labels_support_reads_less_3.2'
            print('False_somatic_DEL.labels_support_reads_less_3.2')
            # return ((dir_name, 'False_somatic_DEL.labels_support_reads_less_3.2'))
        else:
            New_labels = labels[np.where(np.in1d(labels, AimCluster))[0]]
            New_pie = np.array([len(np.where(New_labels == l)[0]) / len(New_labels) for l in AimCluster])
            Newdf['T_or_N'] = Newdf.apply(lambda x: 'N' if 'blood' in x['query_name'] else 'T', axis=1)
            Newdf.to_csv(os.path.join(GaussianMixture_dir, '%s_after_filter.csv' % (name_smple)), sep='\t', index=False)
            csv_filter_path = os.path.join(GaussianMixture_dir, '%s_after_filter.csv' % (name_smple))
            Call_somtic_SV_return = Call_somtic_SV(dir_name,df_all_path)
            if Call_somtic_SV_return[1] == 'Ture_somatic_DEL':
                return_relsult = 'Ture_somatic_DEL'
                print('Ture_somatic_DEL')
                PNG(Newdf,New_x,New_labels,New_pie,means,stds,length_list,labels,png_path)
            else:
                return_relsult = 'False_somatic_DEL.Call_somtic_SV_return'
                # print('False_somatic_DEL.Call_somtic_SV_return[1]')
                # return ((dir_name, 'False_somatic_DEL.Call_somtic_SV_return[1]'))
    return ((dir_name, return_relsult))





### extract DEL loc： bug exist
def get_srm_sv_coord_DEL(chrom, m2r_bam_x, csv_path, region):
    m2r_pysam = pysam.AlignmentFile(m2r_bam_x, 'rb')
    consensus = m2r_pysam.fetch()
    sv_dic = pd.DataFrame(
        columns=['query_name', 'query_sv_start', 'query_sv_end', 'sv_len', 'ref_n', 'sv_ref_start', 'sv_ref_end',
                 'strand'])
    c_x = chrom
    start = int(region.split(':')[1].split('-')[0]) + 500
    end = int(region.split(':')[1].split('-')[1]) - 500
    for AlignedSegmennt in consensus:
        if not AlignedSegmennt.is_secondary:  # 和secondary
            query_n = AlignedSegmennt.query_name
            ref_n = AlignedSegmennt.reference_name
            if ref_n == c_x:  ## Determined to be on that chromosome
                query_sv_start = 1
                query_sv_end = 2
                blocks = AlignedSegmennt.get_blocks()
                for i in range(1, len(blocks)):
                    b = blocks[i]
                    b_last = blocks[i - 1]
                    if (int(b[0]) in range(end - 100, end + 100)) & (int(b_last[1]) in range(start - 100, start + 100)):
                        sv_ref_start = b_last[1]
                        sv_ref_end = b[0]
                        sv_length = int(b[0]) - int(b_last[1])
                        # print(sv_length)
                        if sv_length > 50:
                            sv_dic = sv_dic.append(
                                pd.DataFrame([[query_n, query_sv_start, query_sv_end, sv_length, ref_n,
                                               sv_ref_start,
                                               sv_ref_end, 'positive']],
                                             columns=['query_name', 'query_sv_start', 'query_sv_end',
                                                      'sv_len', 'ref_n', 'sv_ref_start', 'sv_ref_end',
                                                      'strand']), ignore_index=True)
    m2r_pysam.close()
    if sv_dic.shape[0] > 0:
        df = sv_dic[
            ['query_name', 'query_sv_start', 'query_sv_end', 'sv_len', 'ref_n', 'sv_ref_start', 'sv_ref_end', 'strand']]
        # start = int(region.split(':')[1].split('-')[0] ) + 500
        # end = int(region.split(':')[1].split('-')[1]) - 500
        length = int(end) - int(start)
        df['new_length'] = df.apply(lambda x: abs(int(x['sv_len']) - length), axis=1)
        df['new_ref_start'] = df.apply(lambda x: abs(int(x['sv_ref_start']) - int(start)), axis=1)
        zhuangzhi_df = df[(df['new_length'] < 100)]  # & (df['new_ref_start'] < 500) #[(df['new_length'] < 500)]
        zhuangzhi_df.to_csv(csv_path, sep='\t', index=True, header=True)
    else:
        zhuangzhi_df = pd.DataFrame(
            columns=['query_name', 'query_sv_start', 'query_sv_end', 'sv_len', 'ref_n', 'sv_ref_start', 'sv_ref_end',
                     'strand', 'new_length', 'new_ref_start'])
        zhuangzhi_df.to_csv(csv_path, sep='\t', index=True, header=True)




def get_tumor_support_reads_bam(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
    work_dir
    tumor_bam_path
        blood_bam_path
        vcf_region_df: vcf_region_df
        index:
        # wok_dir: bamFile dir
        # SVinfo: SVinfo df
        # index: SV PID
    output:
       {sampel_name_x}_tumor.bam  {sampel_name_x}_tumor_support_reads.fq
    '''
    # region, tumor_bam_path, wok_dir, line_x
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    region = line_x['region']
    chrom = line_x['chr']
    print(dir_name)
    dir_path = os.path.join(wok_dir, dir_name)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    bam_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_tumor.bam'.format(sampel_name_x=dir_name))
    os.system(
        "{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools, bam_ont=tumor_bam_path,
                                                                            region=region, bam_out=bam_out))
    if os.path.getsize(bam_out) == 0:
        os.system("{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools,
                                                                                      bam_ont=blood_bam_path,
                                                                                      region=region, bam_out=bam_out))
    os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    if os.path.exists(bam_out + '.bai'):
        os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    # os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    bed_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_tumor.csv'.format(sampel_name_x=dir_name))
    get_srm_sv_coord_DEL(chrom, bam_out, bed_out, region)


def get_blood_bam(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
        work_dir
        tumor_bam_path
        blood_bam_path
        vcf_region_df: vcf_region_df
        index:
    output:
        record {sampel_name_x}_blood.bam
    '''
    # region, blood_bam_path, wok_dir, dir_name
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    region = line_x['region']
    chrom = line_x['chr']
    dir_path = os.path.join(wok_dir, dir_name)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    bam_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_blood.bam'.format(sampel_name_x=dir_name))
    os.system(
        "{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools, bam_ont=blood_bam_path,
                                                                            region=region, bam_out=bam_out))
    if os.path.getsize(bam_out) == 0:
        os.system("{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools,
                                                                                      bam_ont=blood_bam_path,
                                                                                      region=region, bam_out=bam_out))
    os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    if os.path.exists(bam_out + '.bai'):
        os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    bed_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_blood.csv'.format(sampel_name_x=dir_name))
    get_srm_sv_coord_DEL(chrom, bam_out, bed_out, region)
    # os.system("{bedtools} bamtobed -i {bam_out} -split > {bed_out}".format(bedtools=bedtools, bam_out=bam_out, bed_out=bed_out))


def SelectK(x):
    # choose the best
    print('SelectK')
    print(x)
    if len(x) < 10:
        min_components = 1
        max_components = len(x)
    else:
        min_components = 1
        max_components = 10
    n_components = np.arange(min_components, max_components)  #  Range of groups considered
    models = [GaussianMixture(n, covariance_type='full', random_state=0).fit(x) for n in n_components]
    bics = [m.bic(x) for m in models]  #Calculate the BIC value for each model
    # best_k = n_components[np.argmin(bics)]
    best_gmm = models[np.argmin(bics)]  # Optimal GMM modeling
    labels = best_gmm.predict(x)  # clustering prediction
    return (best_gmm, labels)


# 20231221
def Call_somtic_SV(dir_name, csv_filter_path):
    if os.path.exists(csv_filter_path):
        csv_filter_df = pd.read_csv(csv_filter_path, sep='\t')
        label_list = csv_filter_df['label'].tolist()
        T_label_list = []
        N_label_list = []
        for label in label_list:
            csv_filter_df_label = csv_filter_df[csv_filter_df['label'] == label]
            T_N_list = csv_filter_df_label['T_or_N'].tolist()
            if ('N' not in T_N_list) and (len(T_N_list) >= 3):
                T_label_list.append(label)
            if ('N' in T_N_list):
                N_label_list.append(label)
        csv_filter_df_tumor = csv_filter_df[csv_filter_df['label'].isin(T_label_list)]
        t_test_p = []
        if N_label_list != []:
            # for N_label_list_x in set(N_label_list):
            csv_filter_df_normal_x = csv_filter_df[(csv_filter_df['label'].isin(N_label_list))]
            if (csv_filter_df_normal_x.shape[0] == 1) or (csv_filter_df_tumor.shape[0] == 1):
                normal_sv_len = csv_filter_df_normal_x['sv_len'].tolist()
                normal_sv_len_mean = np.mean(normal_sv_len)
                tumor_sv_len = csv_filter_df_tumor['sv_len'].tolist()
                tumor_sv_len_mean = np.mean(tumor_sv_len)
                if abs(normal_sv_len_mean - tumor_sv_len_mean) < 20:
                    t_test_p = ['False']
            else:
                normal_sv_len = csv_filter_df_normal_x['sv_len'].tolist()
                normal_sv_len_mean = np.mean(normal_sv_len)
                tumor_sv_len = csv_filter_df_tumor['sv_len'].tolist()
                tumor_sv_len_mean = np.mean(tumor_sv_len)
                if abs(normal_sv_len_mean - tumor_sv_len_mean) > 20:
                    p = \
                    stats.ttest_ind(csv_filter_df_normal_x['sv_len'].tolist(), csv_filter_df_tumor['sv_len'].tolist())[
                        1]
                    t_test_p.append(p)
                else:
                    t_test_p = ['False']
            if (csv_filter_df_tumor.shape[0] > 0) and (t_test_p == []):
                csv_filter_df_tumor.to_csv(csv_filter_path.replace('.csv', '.Ture_somatic_SV_DEL.csv'), sep='\t',
                                           index=False, header=True)
                return ((dir_name, 'Ture_somatic_DEL'))
            elif ('False' in t_test_p):
                return ((dir_name, 'False_somatic_DEL.tumor_support_reads_less_3.3'))
            elif ((csv_filter_df_tumor.shape[0] > 0) and (max(t_test_p) < 0.01)):
                csv_filter_df_tumor.to_csv(csv_filter_path.replace('.csv', '.p_value.Ture_somatic_SV_DEL.csv'),
                                           sep='\t', index=False, header=True)
                return ((dir_name, 'Ture_somatic_DEL'))
            else:
                return ((dir_name, 'False_somatic_DEL.tumor_support_reads_less_3.3'))
        elif (N_label_list == []) and csv_filter_df_tumor.shape[0] > 0:
            csv_filter_df_tumor.to_csv(csv_filter_path.replace('.csv', '.Ture_somatic_SV_DEL.csv'), sep='\t',
                                       index=False, header=True)
            return ((dir_name, 'Ture_somatic_DEL'))
        else:
            return ((dir_name, 'False_somatic_DEL.tumor_support_reads_less_3.4'))
    else:
        return ((dir_name, 'False_somatic_DEL.tumor_support_reads_less_3.4'))



###############################################################################################################################
###########  Screening for DEL using Gaussian distribution　　#################
#############################################################################################################################

def get_merge_bam_path(sample_name):
    tumor_bam = os.path.join(os.path.join(somatic_dir, 'tumor', sample_name),
                                     '{DRR258593}_tumor_minimap2.2.26_sorted_tag.bam'.format(DRR258593=sample_name))
    blood_bam = os.path.join(os.path.join(somatic_dir, 'blood', sample_name),
                                     '{DRR258594}_blood_minimap2.2.26_sorted_tag.bam'.format(DRR258594=sample_name))
    return (tumor_bam,blood_bam)



def make_new_df_from_csv(csv_path, sample_name):
    csv_df = pd.read_csv(csv_path, sep='\t')
    csv_df.columns = ['chr', 'start', 'end', 'svid','svtype', 'svlen', 'support_number','number', 'PRECISE', 'STD']
    csv_df_DEL = csv_df.loc[csv_df['svtype'] == 'DEL']
    csv_df_DEL_region = pd.DataFrame()
    csv_df_DEL_region['chr'] = csv_df_DEL['chr']
    csv_df_DEL_region['start'] = csv_df_DEL['start'].astype(int) - 500
    csv_df_DEL_region['end'] = csv_df_DEL['end'].astype(int) + 500
    csv_df_DEL_region['sv_id'] = csv_df_DEL['svid']
    csv_df_DEL_region['sample_name'] = sample_name
    csv_df_DEL_region['region'] = csv_df_DEL_region['chr'].astype(str) + ':' + csv_df_DEL_region['start'].astype(str) + '-' + csv_df_DEL_region['end'].astype(str)
    # csv_df_DEL_region['support_reads_name'] = csv_df_DEL['svreads']
    csv_df_DEL_region['sv_type'] = csv_df_DEL['svtype']
    csv_df_DEL_region['sv_length'] = csv_df_DEL['svlen']
    csv_df_DEL_region['dir_name'] = csv_df_DEL_region['sample_name'] + '_' + csv_df_DEL_region['sv_id'] + '_' + csv_df_DEL_region['chr'] + '_' + csv_df_DEL_region['start'].astype(str) + '_'+csv_df_DEL_region['sv_length'].astype(str)
    csv_df_DEL_region.index = np.array(csv_df_DEL_region['dir_name']) #设置index为dir_name
    return csv_df_DEL_region





