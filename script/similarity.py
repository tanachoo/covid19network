# Author: yoshi
# Date: 6/16/2020
# Project: covid19
# Script: measure similarity for ECv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def main():
    # read biopsy biggest component sif file
    data = pd.read_table("../ECv_result/biopsy/deltaECv_table_Series15_th2.3_biggest.sif", sep='\t', names=('node1','edgetype','node2'))
    edges = []
    for k, l in zip(data['node1'],data['node2']):
        pair = k + '_' + l
        edges.append(pair)

    ecv = pd.read_table('../ECv_result/biopsy/ECv_GSE147507_Log2RPMNormCounts_biopsy.txt', sep='\t', index_col=2)
    ecv_big = ecv.loc[edges,:]
    ecv_big['Series15_COVID19Lung_mean'] = (ecv_big['Series15_COVID19Lung_2']+ecv_big['Series15_COVID19Lung_1'])/2
    ecv_big_ = ecv_big.loc[:,['Series15_HealthyLungBiopsy_2','Series15_HealthyLungBiopsy_1','Series15_COVID19Lung_mean']]
    ecv_big_ = ecv_big_.rename(columns={'Series15_HealthyLungBiopsy_2': 'S15_HealthyLungBiopsy_2',
                                        'Series15_HealthyLungBiopsy_1': 'S15_HealthyLungBiopsy_1',
                                        'Series15_COVID19Lung_mean': 'S15_COVID19Lung'})
    # res = ecv_big_.corr()
    # plt.scatter(ecv_big_['Series15_HealthyLungBiopsy_2'],ecv_big_['Series15_HealthyLungBiopsy_1'])
    # plt.scatter(ecv_big_['Series15_HealthyLungBiopsy_2'],ecv_big_['Series15_COVID19Lung_mean'])


    ecv_vitro = pd.read_table("../ECv_result/vitroRNAseq/network_th005/ECv_network_result_500k_0.05_processing_ParentChild.txt", sep='\t', index_col=0)

    # pick shared-Parent-Child index
    ecv_integ = ecv_big_.join(ecv_vitro, how='inner')

    # res = ecv_integ.corr()
    # sns.heatmap(res, square=True, vmax=1, vmin=-1, center=0)
    # plt.xticks(fontsize=5)
    # plt.yticks(fontsize=5)

    ecv_integ['S1_mock'] = (ecv_integ['Series1_NHBE_Mock_1']
                          + ecv_integ['Series1_NHBE_Mock_2']
                          + ecv_integ['Series1_NHBE_Mock_3'])/3

    ecv_integ['S1_SARS-CoV-2-infected'] = (ecv_integ['Series1_NHBE_SARS-CoV-2_1']
                                         + ecv_integ['Series1_NHBE_SARS-CoV-2_2']
                                         + ecv_integ['Series1_NHBE_SARS-CoV-2_3'])/3

    ecv_integ['S2_mock'] = (ecv_integ['Series2_A549_Mock_1']
                          + ecv_integ['Series2_A549_Mock_2']
                          + ecv_integ['Series2_A549_Mock_3'])/3

    ecv_integ['S2_SARS-CoV-2-infected'] = (ecv_integ['Series2_A549_SARS-CoV-2_1'] 
                                         + ecv_integ['Series2_A549_SARS-CoV-2_2']
                                         + ecv_integ['Series2_A549_SARS-CoV-2_3'])/3

    ecv_integ['S4_mock'] = (ecv_integ['Series4_A549_Mock_1']
                          + ecv_integ['Series4_A549_Mock_2'])/2
    
    ecv_integ['S4_IAV-infected'] = (ecv_integ['Series4_A549_IAV_1']
                                  + ecv_integ['Series4_A549_IAV_2'])/2

    ecv_integ['S5_mock'] = (ecv_integ['Series5_A549_Mock_1']
                          + ecv_integ['Series5_A549_Mock_2']
                          + ecv_integ['Series5_A549_Mock_3'])/3

    ecv_integ['S5_SARS-CoV-2-infected'] = (ecv_integ['Series5_A549_SARS-CoV-2_1']
                                         + ecv_integ['Series5_A549_SARS-CoV-2_2']
                                         + ecv_integ['Series5_A549_SARS-CoV-2_3'])/3

    ecv_integ['S7_mock'] = (ecv_integ['Series7_Calu3_Mock_1']
                          + ecv_integ['Series7_Calu3_Mock_2']
                          + ecv_integ['Series7_Calu3_Mock_3'])/3

    ecv_integ['S7_SARS-CoV-2-infected'] = (ecv_integ['Series7_Calu3_SARS-CoV-2_1']
                                         + ecv_integ['Series7_Calu3_SARS-CoV-2_2']
                                         + ecv_integ['Series7_Calu3_SARS-CoV-2_3'])/3

    ecv_integ['S8_mock'] = (ecv_integ['Series8_A549_Mock_1']
                          + ecv_integ['Series8_A549_Mock_2']
                          + ecv_integ['Series8_A549_Mock_3'])/3

    
    ecv_integ['S8_RSV-infected'] = (ecv_integ['Series8_A549_RSV_1']
                                  + ecv_integ['Series8_A549_RSV_2']
                                  + ecv_integ['Series8_A549_RSV_3'])/3


    ecv_integ['S8_HPIV3-infected'] = (ecv_integ['Series8_A549_HPIV3_1']
                                    + ecv_integ['Series8_A549_HPIV3_2']
                                    + ecv_integ['Series8_A549_HPIV3_3'])/3

    '''
    ecv_integ_ = ecv_integ.loc[:,['S15_HealthyLungBiopsy_2',
                                  'S15_HealthyLungBiopsy_1',
                                  'S15_COVID19Lung',
                                  'S1_mock',
                                  'S1_SARS-CoV-2-infected',
                                  'S2_mock',
                                  'S2_SARS-CoV-2-infected',
                                  'S5_mock',
                                  'S5_SARS-CoV-2-infected',
                                  'S7_mock',
                                  'S7_SARS-CoV-2-infected',
                                  'S4_mock',
                                  'S4_IAV-infected',
                                  'S8_mock',
                                  'S8_RSV-infected',
                                  'S8_HPIV3-infected']]
    '''

    # order: NHBE SARS, Calu3 SARS, A549 MOI0.2 SARS, A549 MOI2 SARS, A549 IAV, A549 RSV, A549 HPIV3
    ecv_integ_ = ecv_integ.loc[:,['S15_COVID19Lung',
                                  'S1_SARS-CoV-2-infected',
                                  'S7_SARS-CoV-2-infected',
                                  'S5_SARS-CoV-2-infected',
                                  'S2_SARS-CoV-2-infected',
                                  'S4_IAV-infected',
                                  'S8_RSV-infected',
                                  'S8_HPIV3-infected']]
    print(ecv_integ_.columns)
    # calculate dist
    #samplenames = ecv_integ_.columns
    samplenames = ['COVID-19 \n (biopsy)',
                   'SARS-CoV-2 \n (NHBE)',
                   'SARS-CoV-2 \n (Calu-3)',
                   'SARS-CoV-2 \n (A549, high)',
                   'SARS-CoV-2 \n (A549, low)',
                   'IAV \n (A549)',
                   'RSV \n (A549)',
                   'HPIV3 \n (A549)']
    dist_matrix = squareform(pdist(ecv_integ_.T.values, metric='cosine')) #metrics: euclidean, cosine, correlation
    #print(dist_matrix)
    df = pd.DataFrame(dist_matrix, index=samplenames, columns=samplenames)
    plt.figure(figsize=(8, 7))
    sns.heatmap(df, square=True, cmap='cividis_r', cbar_kws={'label': 'similarity'})
    plt.xticks(fontsize=12, weight='bold',rotation=70)
    plt.yticks(fontsize=12, weight='bold')
    plt.tight_layout()
    plt.savefig('../ECv_result/biopsy/dist_heatmap_viruses_cosine_TEST.png', dpi=300, format='png')
    plt.close()

    embed_input = ecv_integ_.T.values
    ## Standard scaling for embedding input, mean:0, variance:1
    #scaler = StandardScaler()
    #scaler.fit(embed_input)
    #embed_input = scaler.transform(embed_input)

    '''
    # tSNE
    n_components = 2
    perplexity = 50
    tsne = TSNE(n_components=n_components, init='random', random_state=1111, perplexity=perplexity)
    tsne_output = tsne.fit_transform(embed_input)
    tsne_embed = pd.DataFrame(tsne_output, columns=['tSNE1', 'tSNE2'], index=samplenames)
    plt.scatter(tsne_embed['tSNE1'], tsne_embed['tSNE2'], alpha=1, c='magenta', s=30)
    for x, y, name in zip(tsne_embed['tSNE1'], tsne_embed['tSNE2'], samplenames):
        plt.text(x, y, name, alpha=1, fontsize=10, horizontalalignment='center')
    plt.title('tSNE plot')
    plt.xlabel('tSNE1')
    plt.ylabel('tSNE2')
    plt.savefig('../ECv_result/biopsy/tsne.png', dpi=300, format='png')
    plt.clf()

    # PCA
    n_components = 2
    pca = PCA(n_components)
    pca.fit(embed_input)
    pca_output = pca.transform(embed_input)
    pca_embed = pd.DataFrame(pca_output, columns=['PC1', 'PC2'], index=samplenames)
    plt.scatter(pca_embed['PC1'], pca_embed['PC2'], alpha=1, c='blue', s=30)
    for x, y, name in zip(pca_embed['PC1'], pca_embed['PC2'], samplenames):
        plt.text(x, y, name, alpha=1, fontsize=10, horizontalalignment='center')
    plt.title('PCA plot')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig('../ECv_result/biopsy/pca.png', dpi=300, format='png')
    plt.clf()
    '''

if __name__ == '__main__':
    main()

