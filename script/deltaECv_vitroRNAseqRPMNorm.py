# Auth
# Date: 5/28/2020
# Updated: 5/28/2020
# Project: COVID19
# Script: delta ECv calculation for vitroRNAseq

import pandas as pd
import matplotlib.pyplot as plt

def main():
    filename = '../ECv_result/vitroRNAseqRPMNorm/network_th01/ECv_GSE147507_vitroRNAseqRPMNorm_T500K_1_th010_processing.txt'
    print(f'[LOAD]: {filename}')
    ecv = pd.read_table(filename, sep='\t', header=0)

    ## Calculate deltaECv
    # Series1
    ecv['deltaECv_Series1'] = abs((ecv['Series1_NHBE_Mock_1']
                                 + ecv['Series1_NHBE_Mock_2']
                                 + ecv['Series1_NHBE_Mock_3'])/3
                                - (ecv['Series1_NHBE_SARS-CoV-2_1']
                                 + ecv['Series1_NHBE_SARS-CoV-2_2']
                                 + ecv['Series1_NHBE_SARS-CoV-2_3'])/3)
    # Series2
    ecv['deltaECv_Series2'] = abs((ecv['Series2_A549_Mock_1']
                                 + ecv['Series2_A549_Mock_2']
                                 + ecv['Series2_A549_Mock_3'])/3
                                - (ecv['Series2_A549_SARS-CoV-2_1']
                                 + ecv['Series2_A549_SARS-CoV-2_2']
                                 + ecv['Series2_A549_SARS-CoV-2_3'])/3)
    # Series3
    ecv['deltaECv_Series3'] = abs((ecv['Series3_A549_Mock_1']
                                 + ecv['Series3_A549_Mock_2'])/2
                                - (ecv['Series3_A549_RSV_1']
                                 + ecv['Series3_A549_RSV_2'])/2)
    # Series4
    ecv['deltaECv_Series4'] = abs((ecv['Series4_A549_Mock_1']
                                 + ecv['Series4_A549_Mock_2'])/2
                                - (ecv['Series4_A549_IAV_1']
                                 + ecv['Series4_A549_IAV_2'])/2)
    # Series5
    ecv['deltaECv_Series5'] = abs((ecv['Series5_A549_Mock_1']
                                 + ecv['Series5_A549_Mock_2']
                                 + ecv['Series5_A549_Mock_3'])/3
                                - (ecv['Series5_A549_SARS-CoV-2_1']
                                 + ecv['Series5_A549_SARS-CoV-2_2']
                                 + ecv['Series5_A549_SARS-CoV-2_3'])/3)
    # Series6
    ecv['deltaECv_Series6'] = abs((ecv['Series6_A549-ACE2_Mock_1']
                                 + ecv['Series6_A549-ACE2_Mock_2']
                                 + ecv['Series6_A549-ACE2_Mock_3'])/3
                                - (ecv['Series6_A549-ACE2_SARS-CoV-2_1']
                                 + ecv['Series6_A549-ACE2_SARS-CoV-2_2']
                                 + ecv['Series6_A549-ACE2_SARS-CoV-2_3'])/3)
    # Series7
    ecv['deltaECv_Series7'] = abs((ecv['Series7_Calu3_Mock_1']
                                 + ecv['Series7_Calu3_Mock_2']
                                 + ecv['Series7_Calu3_Mock_3'])/3
                                - (ecv['Series7_Calu3_SARS-CoV-2_1']
                                 + ecv['Series7_Calu3_SARS-CoV-2_1']
                                 + ecv['Series7_Calu3_SARS-CoV-2_1'])/3)
    # Series8
    ecv['deltaECv_Series8_1'] = abs((ecv['Series8_A549_Mock_1']
                                   + ecv['Series8_A549_Mock_2']
                                   + ecv['Series8_A549_Mock_3'])/3
                                  - (ecv['Series8_A549_RSV_1']
                                   + ecv['Series8_A549_RSV_2']
                                   + ecv['Series8_A549_RSV_3'])/3)

    ecv['deltaECv_Series8_2'] = abs((ecv['Series8_A549_Mock_1']
                                   + ecv['Series8_A549_Mock_2']
                                   + ecv['Series8_A549_Mock_3'])/3
                                  - (ecv['Series8_A549_HPIV3_1']
                                   + ecv['Series8_A549_HPIV3_2']
                                   + ecv['Series8_A549_HPIV3_3'])/3)

    # Series9
    ecv['deltaECv_Series9_1'] = abs((ecv['Series9_NHBE_Mock_1']
                                   + ecv['Series9_NHBE_Mock_2']
                                   + ecv['Series9_NHBE_Mock_3']
                                   + ecv['Series9_NHBE_Mock_4'])/4
                                  - (ecv['Series9_NHBE_IAV_1']
                                   + ecv['Series9_NHBE_IAV_2']
                                   + ecv['Series9_NHBE_IAV_3']
                                   + ecv['Series9_NHBE_IAV_4'])/4)

    ecv['deltaECv_Series9_2'] = abs((ecv['Series9_NHBE_Mock_1']
                                   + ecv['Series9_NHBE_Mock_2']
                                   + ecv['Series9_NHBE_Mock_3']
                                   + ecv['Series9_NHBE_Mock_4'])/4
                                  - (ecv['Series9_NHBE_IAVdNS1_1']
                                   + ecv['Series9_NHBE_IAVdNS1_2']
                                   + ecv['Series9_NHBE_IAVdNS1_3']
                                   + ecv['Series9_NHBE_IAVdNS1_4'])/4)

    ecv['deltaECv_Series9_3'] = abs((ecv['Series9_NHBE_Mock_1']
                                   + ecv['Series9_NHBE_Mock_2']
                                   + ecv['Series9_NHBE_Mock_3']
                                   + ecv['Series9_NHBE_Mock_4'])/4
                                  - (ecv['Series9_NHBE_IFNB_4h_1']
                                   + ecv['Series9_NHBE_IFNB_4h_2'])/2)

    ecv['deltaECv_Series9_4'] = abs((ecv['Series9_NHBE_Mock_1']
                                   + ecv['Series9_NHBE_Mock_2']
                                   + ecv['Series9_NHBE_Mock_3']
                                   + ecv['Series9_NHBE_Mock_4'])/4
                                  - (ecv['Series9_NHBE_IFNB_6h_1']
                                   + ecv['Series9_NHBE_IFNB_6h_2'])/2)

    ecv['deltaECv_Series9_5'] = abs((ecv['Series9_NHBE_Mock_1']
                                   + ecv['Series9_NHBE_Mock_2']
                                   + ecv['Series9_NHBE_Mock_3']
                                   + ecv['Series9_NHBE_Mock_4'])/4
                                  - (ecv['Series9_NHBE_IFNB_12h_1']
                                   + ecv['Series9_NHBE_IFNB_12h_2'])/2)

    # Series16
    ecv['deltaECv_Series16_1'] = abs((ecv['Series16_A549-ACE2_Mock_1']
                                    + ecv['Series16_A549-ACE2_Mock_2']
                                    + ecv['Series16_A549-ACE2_Mock_3'])/3
                                   - (ecv['Series16_A549-ACE2_SARS-CoV-2_1']
                                    + ecv['Series16_A549-ACE2_SARS-CoV-2_2']
                                    + ecv['Series16_A549-ACE2_SARS-CoV-2_3'])/3)

    ecv['deltaECv_Series16_2'] = abs((ecv['Series16_A549-ACE2_Mock_1']
                                    + ecv['Series16_A549-ACE2_Mock_2']
                                    + ecv['Series16_A549-ACE2_Mock_3'])/3
                                   - (ecv['Series16_A549-ACE2_SARS-CoV-2_Rux_1']
                                    + ecv['Series16_A549-ACE2_SARS-CoV-2_Rux_2']
                                    + ecv['Series16_A549-ACE2_SARS-CoV-2_Rux_3'])/3)

    ## delta ECv distribution
    delta_Series1 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series1']]
    plt.hist(delta_Series1['deltaECv_Series1'], log=True, color='rosybrown', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series1.png', dpi=300, format='png')
    plt.clf()

    delta_Series2 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series2']]
    plt.hist(delta_Series2['deltaECv_Series2'], log=True, color='darksalmon', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series2.png', dpi=300, format='png')
    plt.clf()

    delta_Series3 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series3']]
    plt.hist(delta_Series3['deltaECv_Series3'], log=True, color='sandybrown', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series3.png', dpi=300, format='png')
    plt.clf()

    delta_Series4 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series4']]
    plt.hist(delta_Series4['deltaECv_Series4'], log=True, color='darkkhaki', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series4.png', dpi=300, format='png')
    plt.clf()

    delta_Series5 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series5']]
    plt.hist(delta_Series5['deltaECv_Series5'], log=True, color='olivedrab', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series5.png', dpi=300, format='png')
    plt.clf()

    delta_Series6 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series6']]
    plt.hist(delta_Series6['deltaECv_Series6'], log=True, color='chartreuse', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series6.png', dpi=300, format='png')
    plt.clf()

    delta_Series7 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series7']]
    plt.hist(delta_Series7['deltaECv_Series7'], log=True, color='seagreen', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series7.png', dpi=300, format='png')
    plt.clf()

    delta_Series8_1 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series8_1']]
    plt.hist(delta_Series8_1['deltaECv_Series8_1'], log=True, color='lightseagreen', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series8_1.png', dpi=300, format='png')
    plt.clf()

    delta_Series8_2 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series8_2']]
    plt.hist(delta_Series8_2['deltaECv_Series8_2'], log=True, color='darkcyan', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series8_2.png', dpi=300, format='png')
    plt.clf()

    delta_Series9_1 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series9_1']]
    plt.hist(delta_Series9_1['deltaECv_Series9_1'], log=True, color='deepskyblue', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series9_1.png', dpi=300, format='png')
    plt.clf()

    delta_Series9_2 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series9_2']]
    plt.hist(delta_Series9_2['deltaECv_Series9_2'], log=True, color='slategray', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series9_2.png', dpi=300, format='png')
    plt.clf()

    delta_Series9_3 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series9_3']]
    plt.hist(delta_Series9_3['deltaECv_Series9_3'], log=True, color='royalblue', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series9_3.png', dpi=300, format='png')
    plt.clf()

    delta_Series9_4 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series9_4']]
    plt.hist(delta_Series9_4['deltaECv_Series9_4'], log=True, color='darkturquoise', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series9_4.png', dpi=300, format='png')
    plt.clf()

    delta_Series9_5 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series9_5']]
    plt.hist(delta_Series9_5['deltaECv_Series9_5'], log=True, color='mediumpurple', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series9_5.png', dpi=300, format='png')
    plt.clf()

    delta_Series16_1 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series16_1']]
    plt.hist(delta_Series16_1['deltaECv_Series16_1'], log=True, color='mediumvioletred', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series16_1.png', dpi=300, format='png')
    plt.clf()

    delta_Series16_2 = ecv.loc[:,['Parent', 'Child', 'deltaECv_Series16_2']]
    plt.hist(delta_Series16_2['deltaECv_Series16_2'], log=True, color='palevioletred', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    plt.savefig('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_histogram_Series16_2.png', dpi=300, format='png')
    plt.clf()


    ## delta ECV threshold
    deltaECv_th = 1

    delta_Series1_th = delta_Series1[delta_Series1['deltaECv_Series1'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series1_th1.txt', 'w') as f:
        delta_Series1_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series2_th = delta_Series2[delta_Series2['deltaECv_Series2'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series2_th1.txt', 'w') as f:
        delta_Series2_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series3_th = delta_Series3[delta_Series3['deltaECv_Series3'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series3_th1.txt', 'w') as f:
        delta_Series3_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series4_th = delta_Series4[delta_Series4['deltaECv_Series4'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series4_th1.txt', 'w') as f:
        delta_Series4_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series5_th = delta_Series5[delta_Series5['deltaECv_Series5'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series5_th1.txt', 'w') as f:
        delta_Series5_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series6_th = delta_Series6[delta_Series6['deltaECv_Series6'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series6_th1.txt', 'w') as f:
        delta_Series6_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series7_th = delta_Series7[delta_Series7['deltaECv_Series7'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series7_th1.txt', 'w') as f:
        delta_Series7_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series8_1_th = delta_Series8_1[delta_Series8_1['deltaECv_Series8_1'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series8_1_th1.txt', 'w') as f:
        delta_Series8_1_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series8_2_th = delta_Series8_2[delta_Series8_2['deltaECv_Series8_2'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series8_2_th1.txt', 'w') as f:
        delta_Series8_2_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series9_1_th = delta_Series9_1[delta_Series9_1['deltaECv_Series9_1'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series9_1_th1.txt', 'w') as f:
        delta_Series9_1_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series9_2_th = delta_Series9_2[delta_Series9_2['deltaECv_Series9_2'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series9_2_th1.txt', 'w') as f:
        delta_Series9_2_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series9_3_th = delta_Series9_3[delta_Series9_3['deltaECv_Series9_3'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series9_3_th1.txt', 'w') as f:
        delta_Series9_3_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series9_4_th = delta_Series9_4[delta_Series9_4['deltaECv_Series9_4'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series9_4_th1.txt', 'w') as f:
        delta_Series9_4_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series9_5_th = delta_Series9_1[delta_Series9_5['deltaECv_Series9_5'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series9_5_th1.txt', 'w') as f:
        delta_Series9_5_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series16_1_th = delta_Series16_1[delta_Series16_1['deltaECv_Series16_1'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series16_1_th1.txt', 'w') as f:
        delta_Series16_1_th.to_csv(f, sep='\t', header=True, index=False)

    delta_Series16_2_th = delta_Series16_2[delta_Series16_2['deltaECv_Series16_2'] >= deltaECv_th]
    with open('../ECv_result/vitroRNAseqRPMNorm/network_th01/deltaECv_table_Series16_2_th1.txt', 'w') as f:
        delta_Series16_2_th.to_csv(f, sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()

