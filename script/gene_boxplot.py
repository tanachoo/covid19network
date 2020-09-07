# Author: yoshi
# Date: 6/19/2020
# Project: covid19
# Script: To make barplot using RNAseq expression data

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def main():
    ## Load raw data
    data = pd.read_table("../data/Blanco_RNAseqData/raw/GSE147507_RawReadCounts_Human.tsv", sep='\t', header=0, index_col=0)
    data_=data.loc[:,['Series2_A549_Mock_1',
                      'Series2_A549_Mock_2',
                      'Series2_A549_Mock_3',
                      'Series2_A549_SARS-CoV-2_1',
                      'Series2_A549_SARS-CoV-2_2',
                      'Series2_A549_SARS-CoV-2_3',
                      'Series5_A549_Mock_1',
                      'Series5_A549_Mock_2',
                      'Series5_A549_Mock_3',
                      'Series5_A549_SARS-CoV-2_1',
                      'Series5_A549_SARS-CoV-2_2',
                      'Series5_A549_SARS-CoV-2_3']]

    ## RPM normalization
    RPMNorm_data = pd.DataFrame(index=data_.index)
    for i in data_.columns:
        RPMNorm_data[i] = data_[i]*1000000/(data_[i].sum())
        log2RPMNorm_data = np.log2(RPMNorm_data+1)

    ## Pick ACE2 (whatever you want)
    ace2 = log2RPMNorm_data.loc[['ACE2'],:]
    ace2 = ace2.reset_index() # reset index
    ace2_ = ace2.melt(id_vars='index', var_name='sample', value_name='RPM normalized counts')
    ace2_['replicates']=['S2_Mock',
                         'S2_Mock',
                         'S2_Mock',
                         'S2_SARS-CoV-2',
                         'S2_SARS-CoV-2',
                         'S2_SARS-CoV-2',
                         'S5_Mock',
                         'S5_Mock',
                         'S5_Mock',
                         'S5_SARS-CoV-2',
                         'S5_SARS-CoV-2',
                         'S5_SARS-CoV-2']

    sns.barplot(x='replicates', y='RPM normalized counts', data=ace2_, palette='Accent', errwidth=0.8, capsize=0.1)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_barplot_ACE2.png', dpi=300, format='png')
    plt.clf()

if __name__ == '__main__':
    main()

