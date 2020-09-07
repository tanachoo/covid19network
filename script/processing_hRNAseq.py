# Author: yoshi
# Date: 5/1/2020
# Project: COVID19 network
# Script: check data quality and process data for BN input

import matplotlib.pyplot as plt
from pandas import plotting
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def main()
    filename = '../Blanco_RNAseqData/GSE147507_RawReadCounts_Human.tsv'
    data = pd.read_table(filename, sep='\t', header=0, index_col=0)
    data_log = np.log2(data+1) # float64
    # data_log = data_log.astype('float32') # if need
    sample_names = list(data.columns)

    ## Check correlation between biological replicates
    series1_mock = data_log[['Series1_NHBE_Mock_1', 'Series1_NHBE_Mock_2', 'Series1_NHBE_Mock_3']]
    plotting.scatter_matrix(series1_mock, diagonal='hist', alpha=0.2, s=1, c='magenta')

    ## Standard scaling
    scaler = StandardScaler()
    scaler.fit(data_log)
    data_log_std = scaler.transform(data_log)

    ## PCA
    n_components = 2
    pca = PCA(n_components)
    pca.fit(data_log_std)
    pca_output = pca.transform(data_log_std)
    embed = pd.DataFrame(pca_output, columns=['PC1', 'PC2'], index=sample_names)

    plt.scatter(embed['PC1'], embed['PC2'], alpha=0.6, c='magenta', s=10)
    for x, y, name in zip(embed['PC1'], embed['PC2'], sample_names):
        plt.text(x, y, name, alpha=0.6, s=5)


    ## RNAseq data curation for BNinput
    mean_threshold = 0.30 
    mean_percentile_threshold = np.percentile(data_log.mean(axis=1), 100*mean_threshold)
    indexfer = data_log.mean(axis=1) >= mean_percentile_threshold
    rnaseq_data_curated = data_log.loc[indexfer,:]

    outputfile = ""
    with open(outputfile, 'w') as f:
        rnaseq_data_curated.to_csv(f, sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()

