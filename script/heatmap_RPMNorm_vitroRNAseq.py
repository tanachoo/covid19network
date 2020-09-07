# Author: yoshi
# Date: 7/8/2020
# Project: COVID-19 network
# Script: heatmap for RNAseq data (fold change)

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def main():
    filename = '../data/Blanco_RNAseqData/raw/GSE147507_RawReadCounts_Human.tsv'
    print(f'\n[LOAD]: {filename}')
    data = pd.read_table(filename, sep='\t', header=0, index_col=0)
    log2RawCounts_data = np.log2(data+1)
    log2RawCounts_data = log2RawCounts_data.drop(columns=['Series15_HealthyLungBiopsy_2',
                                                           'Series15_HealthyLungBiopsy_1',
                                                           'Series15_COVID19Lung_2',
                                                           'Series15_COVID19Lung_1'])
    mean_threshold = 0.30
    mean_percentile_threshold = np.percentile(log2RawCounts_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RawCounts_data.mean(axis=1) >= mean_percentile_threshold
    log2RawCounts_data_ = log2RawCounts_data.loc[indexfer,:]


    ## RPM normalization
    RPMNorm_data = pd.DataFrame(index=data.index)
    for i in data.columns:
        RPMNorm_data[i] = data[i]*1000000/(data[i].sum())
    log2RPMNorm_data = np.log2(RPMNorm_data+1)
    log2RPMNorm_data = log2RPMNorm_data.drop(columns=['Series15_HealthyLungBiopsy_2',
                                                      'Series15_HealthyLungBiopsy_1',
                                                      'Series15_COVID19Lung_2',
                                                      'Series15_COVID19Lung_1'])
    mean_threshold = 0.40
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    log2RPMNorm_data_ = log2RPMNorm_data.loc[indexfer,:]

    FCdata = pd.DataFrame()

    FCdata['S1'] = (log2RPMNorm_data_['Series1_NHBE_SARS-CoV-2_1'] \
                  + log2RPMNorm_data_['Series1_NHBE_SARS-CoV-2_2'] \
                  + log2RPMNorm_data_['Series1_NHBE_SARS-CoV-2_3'])/3 \
                 - (log2RPMNorm_data_['Series1_NHBE_Mock_1'] \
                  + log2RPMNorm_data_['Series1_NHBE_Mock_2'] \
                  + log2RPMNorm_data_['Series1_NHBE_Mock_3'])/3

    FCdata['S2'] = (log2RPMNorm_data_['Series2_A549_SARS-CoV-2_1'] \
                  + log2RPMNorm_data_['Series2_A549_SARS-CoV-2_2'] \
                  + log2RPMNorm_data_['Series2_A549_SARS-CoV-2_3'])/3 \
                 - (log2RPMNorm_data_['Series2_A549_Mock_1'] \
                  + log2RPMNorm_data_['Series2_A549_Mock_2'] \
                  + log2RPMNorm_data_['Series2_A549_Mock_3'])/3

    FCdata['S3'] = (log2RPMNorm_data_['Series3_A549_RSV_1'] \
                  + log2RPMNorm_data_['Series3_A549_RSV_2'])/2 \
                 - (log2RPMNorm_data_['Series3_A549_Mock_1'] \
                  + log2RPMNorm_data_['Series3_A549_Mock_1'])/2

    FCdata['S4'] = (log2RPMNorm_data_['Series4_A549_IAV_1'] \
                  + log2RPMNorm_data_['Series4_A549_IAV_1'])/2 \
                 - (log2RPMNorm_data_['Series4_A549_Mock_1'] \
                  + log2RPMNorm_data_['Series4_A549_Mock_1'])/2

    FCdata['S5'] = (log2RPMNorm_data_['Series5_A549_SARS-CoV-2_1'] \
                  + log2RPMNorm_data_['Series5_A549_SARS-CoV-2_2'] \
                  + log2RPMNorm_data_['Series5_A549_SARS-CoV-2_3'])/3 \
                 - (log2RPMNorm_data_['Series5_A549_Mock_1'] \
                  + log2RPMNorm_data_['Series5_A549_Mock_2'] \
                  + log2RPMNorm_data_['Series5_A549_Mock_3'])/3

    FCdata['S6'] = (log2RPMNorm_data_['Series6_A549-ACE2_SARS-CoV-2_1'] \
                  + log2RPMNorm_data_['Series6_A549-ACE2_SARS-CoV-2_2'] \
                  + log2RPMNorm_data_['Series6_A549-ACE2_SARS-CoV-2_3'])/3 \
                 - (log2RPMNorm_data_['Series6_A549-ACE2_Mock_1'] \
                  + log2RPMNorm_data_['Series6_A549-ACE2_Mock_2'] \
                  + log2RPMNorm_data_['Series6_A549-ACE2_Mock_3'])/3

    FCdata['S7'] = (log2RPMNorm_data_['Series7_Calu3_SARS-CoV-2_1'] \
                  + log2RPMNorm_data_['Series7_Calu3_SARS-CoV-2_2'] \
                  + log2RPMNorm_data_['Series7_Calu3_SARS-CoV-2_3'])/3 \
                 - (log2RPMNorm_data_['Series7_Calu3_Mock_1'] \
                  + log2RPMNorm_data_['Series7_Calu3_Mock_2'] \
                  + log2RPMNorm_data_['Series7_Calu3_Mock_3'])/3

    FCdata['S8_1'] = (log2RPMNorm_data_['Series8_A549_RSV_1'] \
                    + log2RPMNorm_data_['Series8_A549_RSV_2'] \
                    + log2RPMNorm_data_['Series8_A549_RSV_3'])/3 \
                   - (log2RPMNorm_data_['Series8_A549_Mock_1'] \
                    + log2RPMNorm_data_['Series8_A549_Mock_2'] \
                    + log2RPMNorm_data_['Series8_A549_Mock_3'])/3 \

    FCdata['S8_2'] = (log2RPMNorm_data_['Series8_A549_HPIV3_1'] \
                    + log2RPMNorm_data_['Series8_A549_HPIV3_2'] \
                    + log2RPMNorm_data_['Series8_A549_HPIV3_3'])/3 \
                   - (log2RPMNorm_data_['Series8_A549_Mock_1'] \
                    + log2RPMNorm_data_['Series8_A549_Mock_2'] \
                    + log2RPMNorm_data_['Series8_A549_Mock_3'])/3

    FCdata['S9_1'] = (log2RPMNorm_data_['Series9_NHBE_IAV_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_IAV_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_IAV_3'])/3 \
                   - (log2RPMNorm_data_['Series9_NHBE_Mock_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_3'])/3

    FCdata['S9_2'] = (log2RPMNorm_data_['Series9_NHBE_IAVdNS1_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_IAVdNS1_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_IAVdNS1_3'])/3 \
                   - (log2RPMNorm_data_['Series9_NHBE_Mock_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_3'])/3

    FCdata['S9_3'] = (log2RPMNorm_data_['Series9_NHBE_IFNB_4h_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_IFNB_4h_2'])/2 \
                   - (log2RPMNorm_data_['Series9_NHBE_Mock_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_3'])/3

    FCdata['S9_4'] = (log2RPMNorm_data_['Series9_NHBE_IFNB_6h_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_IFNB_6h_2'])/2 \
                   - (log2RPMNorm_data_['Series9_NHBE_Mock_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_3'])/3

    FCdata['S9_5'] = (log2RPMNorm_data_['Series9_NHBE_IFNB_12h_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_IFNB_12h_2'])/2 \
                   - (log2RPMNorm_data_['Series9_NHBE_Mock_1'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_2'] \
                    + log2RPMNorm_data_['Series9_NHBE_Mock_3'])/3

    FCdata['S16_1'] = (log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_1'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_2'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_3'])/3 \
                    - (log2RPMNorm_data_['Series16_A549-ACE2_Mock_1'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_Mock_2'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_Mock_3'])/3

    FCdata['S16_2'] = (log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_Rux_1'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_Rux_2'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_SARS-CoV-2_Rux_3'])/3 \
                    - (log2RPMNorm_data_['Series16_A549-ACE2_Mock_1'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_Mock_2'] \
                     + log2RPMNorm_data_['Series16_A549-ACE2_Mock_3'])/3

    '''
    print(f'FCdata shape: {FCdata.shape}')   
    sns.clustermap(FCdata, method='ward', metric='euclidean',
                   yticklabels=False, xticklabels=True,
                   row_cluster=True, col_cluster=False,
                   cmap='seismic', vmax=10, vmin=-10, center=0)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCountsFC_vitroRNAseq_heatmap_TEST.png', dpi=300, format='png')
    plt.close()
    '''

    '''
    print(f'log2RPMNorm matrix: {log2RPMNorm_data_.shape}')
    sns.clustermap(log2RPMNorm_data_, method='ward', metric='euclidean', z_score=1,
                   yticklabels=False, xticklabels=False,
                   row_cluster=True, col_cluster=True,
                   cmap='seismic', vmax=6, vmin=-6, center=0)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_vitroRNAseq_heatmap_TEST.png', dpi=300, format='png')
    plt.close()
    '''


    print(f'log2RawCounts matrix: {log2RawCounts_data_.shape}')
    sns.clustermap(log2RawCounts_data_, method='ward', metric='euclidean', z_score=1,
                   yticklabels=False, xticklabels=False,
                   row_cluster=True, col_cluster=True,
                   cmap='PRGn_r', vmax=6, vmin=-6, center=0)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawCounts_vitroRNAseq_heatmap_TEST.png', dpi=300, format='png')
    plt.close()


if __name__ == '__main__':
    main()


