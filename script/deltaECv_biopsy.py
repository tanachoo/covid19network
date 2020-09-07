# Auth
# Date: 6/12/2020
# Updated: 6/12/2020
# Project: COVID19
# Script: delta ECv calculation for biopsy dataset (Series15)

import pandas as pd
import matplotlib.pyplot as plt

def main():
    filename = '../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy.txt'
    print(f'[LOAD]: {filename}')
    ecv = pd.read_table(filename, sep='\t', header=0)

    ## Calculate deltaECv
    # Series1
    ecv['deltaECv_Series15'] = abs((ecv['Series15_HealthyLungBiopsy_2']
                                  + ecv['Series15_HealthyLungBiopsy_1'])/2
                                 - (ecv['Series15_COVID19Lung_2']
                                  + ecv['Series15_COVID19Lung_1'])/2)
 
    ## delta ECv distribution
    delta_Series15 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series15']]
    plt.hist(delta_Series15['deltaECv_Series15'], log=True, color='gold', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    #plt.hist(delta_Series15['deltaECv_Series15'], log=True, color='yellow', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    #plt.title('ΔECv histogram')
    plt.xlabel('ΔECv', fontsize=24)
    plt.ylabel('frequency', fontsize=24)
    plt.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig('../ECv_result/biopsy/deltaECv_histogram_Series15_gold.png', dpi=300, format='png')
    plt.close()

    ## delta ECV threshold
    deltaECv_th = 1
    '''
    delta_Series15_th = delta_Series15[delta_Series15['deltaECv_Series15'] >= deltaECv_th]
    with open('../ECv_result/biopsy/deltaECv_table_Series15_th1.txt', 'w') as f:
        delta_Series15_th.to_csv(f, sep='\t', header=True, index=False)
    '''

if __name__ == '__main__':
    main()

