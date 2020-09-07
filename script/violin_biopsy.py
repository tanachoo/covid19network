# Author: yoshi
# Date: 6/16/2020
# Updated: 6/24/2020
# Project: covid19
# Script: violinplot

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    biopsy_ecv = pd.read_table('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy.txt',sep='\t',index_col=2)
    biopsy_ecv_ = biopsy_ecv.loc[:,['Series15_HealthyLungBiopsy_2',
                                    'Series15_HealthyLungBiopsy_1',
                                    'Series15_COVID19Lung_2',
                                    'Series15_COVID19Lung_1']]

    biopsy_ecv_=biopsy_ecv_.rename(columns={'Series15_HealthyLungBiopsy_2':'Healthy_2',
                                            'Series15_HealthyLungBiopsy_1':'Healthy_1',
                                            'Series15_COVID19Lung_2':'COVID19_2',
                                            'Series15_COVID19Lung_1':'COVID19_1'})
    biopsy_ecv_=biopsy_ecv_.loc[:,['Healthy_1', 'Healthy_2', 'COVID19_1', 'COVID19_2']]

    # ECv violinplot
    sns.violinplot(data=biopsy_ecv_, palette='Accent', inner='box') # Pattern 1
    plt.ylabel('ECv')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_violinplot.png', dpi=300, format='png')
    plt.clf()

    # ECv boxplot
    sns.boxplot(data=biopsy_ecv_, palette='Accent', width=0.5) # Pattern 1
    plt.ylabel('ECv')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_boxplot.png', dpi=300, format='png')
    plt.clf()

    '''
    #aa=biopsy_.reset_index() # reset index
    #aa.set_index('name') # set index
    #aaa=aa.melt(id_vars='name',var_name='sample',value_name='ECv')
    #sns.violinplot(x='sample',y='value',data=aaa,palette='Accent', inner='box') # Pattern 2
    '''



    '''
    ## Raw count
    biopsy_rawcount = pd.read_table("../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy.txt", sep='\t', header=0, index_col=0)
    biopsy_rawcount_ = biopsy_rawcount.rename(columns={'Series15_HealthyLungBiopsy_2':'Healthy_2',
                                                       'Series15_HealthyLungBiopsy_1':'Healthy_1',
                                                       'Series15_COVID19Lung_2':'COVID19_2',
                                                       'Series15_COVID19Lung_1':'COVID19_1'})
    biopsy_rawcount_ = biopsy_rawcount_.loc[:,['Healthy_1', 'Healthy_2', 'COVID19_1', 'COVID19_2']]

    sns.violinplot(data=biopsy_rawcount_, palette='Accent', inner='box') # Pattern 1
    plt.ylabel('Raw counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy_violinplot.png', dpi=300, format='png')
    plt.clf()

    sns.boxplot(data=biopsy_rawcount_, palette='Accent', width=0.5) # Pattern 1
    plt.ylabel('Raw counts')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy_boxplot.png', dpi=300, format='png')
    plt.clf()
    '''





if __name__ == '__main__':
    main()

