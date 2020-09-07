# Author: yoshi
# Date: 5/15/2020
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
    data = data.drop(columns=['Series15_HealthyLungBiopsy_2','Series15_HealthyLungBiopsy_1',
                              'Series15_COVID19Lung_2','Series15_COVID19Lung_1'])
    data_log = np.log2(data+1) # float64
    # data_log = data_log.astype('float32') # if need
    sample_names = list(data.columns)
    print(f'matrix shape: {data_log.shape}')

    ## check expression distributions
    data_log.plot.box()
    plt.savefig('../rawRNAseq_boxplot.png', dpi=300, format='png')
    plt.clf()

    ## Check correlation between biological replicates
    #series1_mock = data_log[['Series1_NHBE_Mock_1', 'Series1_NHBE_Mock_2', 'Series1_NHBE_Mock_3']]
    #plotting.scatter_matrix(series1_mock, diagonal='hist', alpha=0.2, s=1, c='magenta')

    ## Standard scaling for PCA input, mean:0, variance:1
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

    ## RNAseq data curation for BN input
    plt.hist(data_log.mean(axis='columns'), alpha=0.5, bins=20)
    plt.savefig('../raw_mean_distribution.png', dpi=300, format='png')
    
    mean_threshold = 0.15 
    mean_percentile_threshold = np.percentile(data_log.mean(axis=1), 100*mean_threshold)
    indexfer = data_log.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = data_log.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.5, bins=20)
    plt.savefig('../th015_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.20 
    mean_percentile_threshold = np.percentile(data_log.mean(axis=1), 100*mean_threshold)
    indexfer = data_log.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = data_log.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.5, bins=20)
    plt.savefig('../th020_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.25 
    mean_percentile_threshold = np.percentile(data_log.mean(axis=1), 100*mean_threshold)
    indexfer = data_log.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = data_log.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.5, bins=20)
    plt.savefig('../th025_mean_distribution.png', dpi=300, format='png')

    mean_threshold = 0.30 
    mean_percentile_threshold = np.percentile(data_log.mean(axis=1), 100*mean_threshold)
    indexfer = data_log.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = data_log.loc[indexfer,:]
    plt.hist(rnaseq_data_curated.mean(axis='columns'), alpha=0.5, bins=20)
    plt.savefig('../th030_mean_distribution.png', dpi=300, format='png')
    

    outputfile = "../GSE147507_vitroRNAseq.txt"
    print(f'[SAVE]: {outputfile}\n')
    with open(outputfile, 'w') as f:
        rnaseq_data_curated.to_csv(f, sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()

