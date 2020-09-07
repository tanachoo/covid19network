# Author: yoshi
# Date: 7/13/2020
# Project: covid19 network
# Script: depict heatmap for biopsy samples

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():

    filename = '../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_forRC2.txt'
    print(f'[LOAD]: {filename}')
    data = pd.read_table(filename, header=0, sep='\t', index_col=0)
    print(f'matrix shape: {data.shape}')

    data = data.drop(columns=['COVID19Lung_mean'])
 
    sns.clustermap(data, method='ward', metric='euclidean', z_score=1,
                   yticklabels=False, xticklabels=False,
                   row_cluster=True, col_cluster=False,
                   cbar_kws={'label': 'RNAseq Z score'},
                   cmap='seismic', vmax=5, vmin=-5, center=0)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_heatmap_TEST.png', dip=300, format='png')
    plt.close()

if __name__ == '__main__':
    main()
