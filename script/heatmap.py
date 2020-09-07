# Author: yoshi
# Date: 6/17/2020
# Project: covid19
# Script: Depict heatmap with clustering

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():

    data = pd.read_table("../ECv_result/vitroRNAseq/network_th005/dealtaECv_Series2and5_unionedges.txt", sep='\t', index_col=0)

    sns.clustermap(data, method='ward', metric='euclidean',
                   yticklabels=False, xticklabels=False,
                   cmap='Oranges', cbar_kws={'label': 'ΔECv'})

    plt.savefig('../ECv_result/vitroRNAseq/network_th005/heatmap/dealtaECv_Series2and5_union_edgesheatmap.png', dpi=300, format='png')

if __name__ == '__main__':
    main()
