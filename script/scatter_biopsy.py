# Author: yoshi
# Date: 6/24/2020
# Updated: 6/25/2020
# Project: covid19
# Script: write satter plot

import pandas as pd
import matplotlib.pyplot as plt

def main():

    rpmdata = pd.read_table('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy.txt', header=0, sep='\t', index_col=0)
    rawdata = pd.read_table('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy.txt', header=0, sep='\t', index_col=0)
    rpmdata_rmZero = pd.read_table('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_forRC2.txt', header=0, sep='\t', index_col=0)
    ecv = pd.read_table('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy.txt', header=0, sep='\t',index_col=2)

    '''
    ## RPM normalized Counts
    # healthy
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata['Series15_HealthyLungBiopsy_1'], rpmdata['Series15_HealthyLungBiopsy_2'], alpha=0.5, s=1, c='green')
    plt.xlabel('healthy 1')
    plt.ylabel('healthy 2')
    plt.title('RPM normarized counts scatter plot (healthy)')
    plt.xlim(-0.4, 14)
    plt.ylim(-0.4, 14)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy.png', dpi=300, format='png')
    plt.clf()
    # covid-19
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata['Series15_COVID19Lung_1'], rpmdata['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='magenta')
    plt.xlabel('COVID-19 1')
    plt.ylabel('COVID-19 2')
    plt.title('RPM normarized counts scatter plot (COVID-19)')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_scatterplot_covid19.png', dpi=300, format='png')
    plt.clf()


    ## RawCounts
    # healthy
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rawdata['Series15_HealthyLungBiopsy_1'], rawdata['Series15_HealthyLungBiopsy_2'], alpha=0.5, s=1, c='green')
    plt.xlabel('healthy 1')
    plt.ylabel('healthy 2')
    plt.title('Raw counts scatter plot (healthy)')
    plt.xlim(-0.4, 18)
    plt.ylim(-0.4, 18)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy_scatterplot_healthy.png', dpi=300, format='png')
    plt.clf()
    # covid-19
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rawdata['Series15_COVID19Lung_1'], rawdata['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='magenta')
    plt.xlabel('COVID-19 1')
    plt.ylabel('COVID-19 2')
    plt.title('Raw counts scatter plot (COVID-19)')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RawReadCounts_biopsy_scatterplot_covid19.png', dpi=300, format='png')
    plt.clf()


    ## RPMnormalized remove ZeroCounts
    # healthy
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_1'], rpmdata_rmZero['Series15_HealthyLungBiopsy_2'], alpha=0.5, s=1, c='green')
    plt.xlabel('healthy 1')
    plt.ylabel('healthy 2')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot (healthy)')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy.png', dpi=300, format='png')
    plt.clf()
    # covid-19
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_COVID19Lung_1'], rpmdata_rmZero['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='magenta')
    plt.xlabel('COVID-19 1')
    plt.ylabel('COVID-19 2')
    plt.title('RPM normalized counts post remove Zero counts scatter plot (COVID-19)')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_covid19.png', dpi=300, format='png')
    plt.clf()
    '''

    '''
    ## crossing pattern healthy VS covid-19
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_1'], rpmdata_rmZero['Series15_COVID19Lung_1'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 1')
    plt.ylabel('COVID-19 1')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy1covid191.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_1'], rpmdata_rmZero['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 1')
    plt.ylabel('COVID-19 2')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy1covid192.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_2'], rpmdata_rmZero['Series15_COVID19Lung_1'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 2')
    plt.ylabel('COVID-19 1')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy2covid191.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_2'], rpmdata_rmZero['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 2')
    plt.ylabel('COVID-19 2')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy2covid192.png', dpi=300, format='png')
    plt.clf()
    '''

    '''
    ##covid19 mean version
    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_1'], rpmdata_rmZero['COVID19Lung_mean'], alpha=0.5, s=1, c='orange')
    plt.xlabel('healthy 1')
    plt.ylabel('COVID-19 mean')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy1covid19mean.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(rpmdata_rmZero['Series15_HealthyLungBiopsy_2'], rpmdata_rmZero['COVID19Lung_mean'], alpha=0.5, s=1, c='orange')
    plt.xlabel('healthy 2')
    plt.ylabel('COVID-19 mean')
    plt.xlim(-0.4, 16)
    plt.ylim(-0.4, 16)
    plt.title('RPM normalized counts post remove Zero counts scatter plot')
    plt.savefig('../data/Blanco_RNAseqData/processing/GSE147507_Log2RPMNormCounts_biopsy_removeZeroCounts_scatterplot_healthy2covid19mean.png', dpi=300, format='png')
    plt.clf()
    '''



    ## ECv
    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_HealthyLungBiopsy_1'], ecv['Series15_HealthyLungBiopsy_2'], alpha=0.5, s=1, c='green')
    plt.xlabel('healthy 1')
    plt.ylabel('healthy 2')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_COVID19Lung_1'], ecv['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='magenta')
    plt.xlabel('COVID-19 1')
    plt.ylabel('COVID-19 2')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_covid19.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_HealthyLungBiopsy_1'], ecv['Series15_COVID19Lung_1'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 1')
    plt.ylabel('COVID-19 1')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy1covid191.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_HealthyLungBiopsy_1'], ecv['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 1')
    plt.ylabel('COVID-19 2')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy1covid192.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_HealthyLungBiopsy_2'], ecv['Series15_COVID19Lung_1'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 2')
    plt.ylabel('COVID-19 1')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy2covid191.png', dpi=300, format='png')
    plt.clf()

    plt.figure(figsize=(8, 8)) 
    plt.scatter(ecv['Series15_HealthyLungBiopsy_2'], ecv['Series15_COVID19Lung_2'], alpha=0.5, s=1, c='blue')
    plt.xlabel('healthy 2')
    plt.ylabel('COVID-19 2')
    plt.xlim(-7, 17)
    plt.ylim(-7, 17)
    plt.title('ECv scatter plot')
    plt.savefig('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy_scatterplot_healthy2covid192.png', dpi=300, format='png')
    plt.clf()

if __name__ == '__main__':
    main()


