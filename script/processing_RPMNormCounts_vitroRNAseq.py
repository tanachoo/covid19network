# Author: yoshi
# Date: 5/28/2020
# Project: COVID19 network
# Script: check data quality and process data for BN input

import matplotlib.pyplot as plt
from pandas import plotting
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def main():
    filename = '../data/Blanco_RNAseqData/raw/GSE147507_RawReadCounts_Human.tsv'
    print(f'\n[LOAD]: {filename}')
    data = pd.read_table(filename, sep='\t', header=0, index_col=0)

    ## RPM normalization
    RPMNorm_data = pd.DataFrame(index=data.index)
    for i in data.columns:
        RPMNorm_data[i] = data[i]*1000000/(data[i].sum())
    log2RPMNorm_data = np.log2(RPMNorm_data+1)
    '''
    outputfile1 = "../data/Blanco_RNAseqData/processing/GSE147507_RPMNormCounts.tsv"
    outputfile2 = "../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts.tsv"
    print(f'[SAVE] RPM normalized file: {outputfile1}\n'
          f'[SAVE] log2 RPM normalized file: {outputfile2}')
    with open(outputfile1, 'w') as f:
        RPMNorm_data.to_csv(f, sep='\t', header=True, index=True)
    with open(outputfile2, 'w') as f:
        log2RPMNorm_data.to_csv(f, sep='\t', header=True, index=True)
    '''

    log2RPMNorm_data = log2RPMNorm_data.drop(columns=['Series15_HealthyLungBiopsy_2', 
                                                      'Series15_HealthyLungBiopsy_1',
                                                      'Series15_COVID19Lung_2',
                                                      'Series15_COVID19Lung_1'])

    print(f'log2 RPM Normalized matrix shape: {log2RPMNorm_data.shape}')

    ## check expression distributions
    plt.figure(figsize=(32,48))
    log2RPMNorm_data.plot.box()
    plt.xticks(rotation=80)
    plt.tick_params(labelsize=2)
    plt.savefig('../log2RPMNormCounts_boxplot.png', dpi=300, format='png')
    plt.clf()

    ## Check correlation between biological replicates
    #series1_mock = data_log[['Series1_NHBE_Mock_1', 'Series1_NHBE_Mock_2', 'Series1_NHBE_Mock_3']]
    #plotting.scatter_matrix(series1_mock, diagonal='hist', alpha=0.2, s=1, c='magenta')

    '''
    ## Standard scaling for PCA input, mean:0, variance:1
    sample_names = list(data.columns)
    scaler = StandardScaler()
    scaler.fit(data_log)
    data_log_std = scaler.transform(data_log) #ndarray
    #data_log_std_df=pd.DataFrame(data_log_std, columns=sample_names) #dataframe

    ## PCA
    n_components = 2
    pca = PCA(n_components)
    pca.fit(data_log_std.T)
    pca_output = pca.transform(data_log_std.T)
    embed = pd.DataFrame(pca_output, columns=['PC1', 'PC2'], index=sample_names)

    plt.scatter(embed['PC1'], embed['PC2'], alpha=0.6, c='magenta', s=10)
    for x, y, name in zip(embed['PC1'], embed['PC2'], sample_names):
        plt.text(x, y, name, alpha=0.6, fontsize=3)
    plt.savefig('../pca.png', dpi=300, format='png')
    plt.clf()
    '''

    ## RNAseq data curation for BN input
    plt.hist(log2RPMNorm_data.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_raw_mean_distribution.png', dpi=300, format='png')
    
    mean_threshold = 0.15 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th015_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.20 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th020_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.25 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th025_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.30 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th030_mean_distribution.png', dpi=300, format='png')
    
    mean_threshold = 0.35 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th035_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.40 
    mean_percentile_threshold = np.percentile(log2RPMNorm_data.mean(axis=1), 100*mean_threshold)
    indexfer = log2RPMNorm_data.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = log2RPMNorm_data.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.3, bins=40)
    plt.savefig('../RPMNorm_th040_mean_distribution.png', dpi=300, format='png')

    outputfile = "../GSE147507_RPMNorm_vitroRNAseq.txt"
    print(f'[SAVE]: {outputfile}\n'
          f'final matrix: {rnaseq_data_curated.shape}\n')
    with open(outputfile, 'w') as f:
        rnaseq_data_curated.to_csv(f, sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()

