# Author: yoshi
# Date: 6/15/2020
# Updated: 6/17/2020
# Project: covid19
# Script: To remove ZeroCounts data from biopsy dataset


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    data = pd.read_table("../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy.txt",
                         sep='\t', header=0, index_col=0)

    data_=data.rename(columns={'Series15_HealthyLungBiopsy_2':'Healthy_2',
                               'Series15_HealthyLungBiopsy_1':'Healthy_1',
                               'Series15_COVID19Lung_2':'COVID19_2',
                               'Series15_COVID19Lung_1':'COVID19_1'})
    data_=data_.loc[:,['Healthy_1', 'Healthy_2', 'COVID19_1', 'COVID19_2']]


    '''
    sns.violinplot(data=data_, palette='Accent', inner='box')
    plt.ylabel('RPM normarized counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_violinplot.png', dpi=300, format='png')
    plt.clf()
    '''

    sns.boxplot(data=data_, palette='Accent', width=0.5)
    plt.ylabel('RPM normarized counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_boxplot.png', dpi=300, format='png')
    plt.clf()


    # remove genes with at least one ZeroCounts RPMnormalized Reads
    data_rmzero = data_[data_ > 0].dropna(how='any')


    '''
    sns.violinplot(data=data_rmzero, palette='Accent', inner='box')
    plt.ylabel('RPM normarized counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_violinplot.png', dpi=300, format='png')
    plt.clf()
    '''

    sns.boxplot(data=data_rmzero, palette='Accent', width=0.5)
    plt.ylabel('RPM normarized counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_boxplot.png', dpi=300, format='png')
    plt.clf()


    '''
    with open('GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts.txt', 'w') as f:
        data_rmzero.to_csv(f, sep='\t', header=True, index=True)
    '''


if __name__ == '__main__':
    main()

