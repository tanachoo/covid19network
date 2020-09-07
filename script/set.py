# Author: yoshi
# Date: 6/17/2020
# Project: covid19
# Script: calculate set elements of deltaECv edges

import pandas as pd

def main():
    ## Load biopsy edges
    biopsy=pd.read_table("../ECv_result/biopsy/deltaECv_table_Series15_th2.3_biggest.sif",sep='\t',names=('node1','edgetype','node2'))
    biopsy_edges=[]
    for k,l in zip(biopsy['node1'],biopsy['node2']):
        pair=k+'_'+l
        biopsy_edges.append(pair)

    ## Load different viral loads Series(2,5)
    series2=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_table_Series2_th1.txt",sep='\t',header=0,index_col=False)
    series2_edges=[]
    for k,l in zip(series2['Parent'],series2['Child']):
        pair=k+'_'+l
        series2_edges.append(pair)

    series5=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_table_Series5_th1.txt",sep='\t',header=0,index_col=False)
    series5_edges=[]
    for k,l in zip(series5['Parent'],series5['Child']):
        pair=k+'_'+l
        series5_edges.append(pair)

    len(set(series5_edges)&set(series2_edges))
    len(set(biopsy_edges)&set(series2_edges))
    len(set(biopsy_edges)&set(series5_edges))
    len(set(biopsy_edges)&set(series5_edges)&set(series2_edges))


    ## Calculate Venn Series2/5/biopsy
    series1u5=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_edgesshared_Series1_Series5.graph.tsv",sep='\t',names=('node1','node2'),index_col=False)

    series1u5_edges=[]
    for k,l in zip(series1u5['node1'],series1u5['node2']):
        pair=k+'_'+l
        series1u5_edges.append(pair)


    series1u7=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_edgesshared_Series1_Series7.graph.tsv",sep='\t',names=('node1','node2'),index_col=False)

    series1u7_edges=[]
    for k,l in zip(series1u7['node1'],series1u7['node2']):
        pair=k+'_'+l
        series1u7_edges.append(pair)


    series5u7=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_edgesshared_Series5_Series7.graph.tsv",sep='\t',names=('node1','node2'),index_col=False)

    series5u7_edges=[]
    for k,l in zip(series5u7['node1'],series5u7['node2']):
        pair=k+'_'+l
        series5u7_edges.append(pair)


    len(set(series1u7_edges)|set(series1u5_edges)|set(series5u7_edges))
    series_cells=list(set(series1u7_edges)|set(series1u5_edges)|set(series5u7_edges))


    ## Load different viruses (Series4/5/8_1/8_2)
    series4=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_table_Series4_th1.txt",sep='\t',header=0,index_col=False)
    series4_edges=[]
    for k,l in zip(series4['Parent'],series4['Child']):
        pair=k+'_'+l
        series4_edges.append(pair)

    series8_1=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_table_Series8_1_th1.txt",sep='\t',header=0,index_col=False)
    series8_1_edges=[]
    for k,l in zip(series8_1['Parent'],series8_1['Child']):
        pair=k+'_'+l
        series8_1_edges.append(pair)

    series8_2=pd.read_table("../ECv_result/vitroRNAseq/network_th005/deltaECv_table_Series8_2_th1.txt",sep='\t',header=0,index_col=False)
    series8_2_edges=[]
    for k,l in zip(series8_2['Parent'],series8_2['Child']):
        pair=k+'_'+l
        series8_2_edges.append(pair)


    ## Calculate Venn
    A6=list(set(series8_1_edges)&set(series8_2_edges)&set(series4_edges)&set(series5_edges))
    A6A7=list(set(series8_1_edges)&set(series8_2_edges)&set(series4_edges))
    A5A6=list(set(series8_1_edges)&set(series8_2_edges)&set(series5_edges))
    A6A10=list(set(series8_1_edges)&set(series4_edges)&set(series5_edges))
    A2A6=list(set(series8_2_edges)&set(series4_edges)&set(series5_edges))

    A5A6A7A8=list(set(series8_2_edges)&set(series8_1_edges))
    A5A6A9A10=list(set(series5_edges)&set(series8_1_edges))
    A2A6A10A14=list(set(series5_edges)&set(series4_edges))
    A6A7A10A11=list(set(series4_edges)&set(series8_1_edges))
    A1A2A5A6=list(set(series5_edges)&set(series8_2_edges))
    A2A3A6A7=list(set(series4_edges)&set(series8_2_edges))

    A7=list(set(A6A7)-set(A6))
    A5=list(set(A5A6)-set(A6))
    A10=list(set(A6A10)-set(A6))
    A2=list(set(A2A6)-set(A6))

    A8=list(set(A5A6A7A8)-set(A6)-set(A5)-set(A7))
    A9=list(set(A5A6A9A10)-set(A6)-set(A5)-set(A10))
    A14=list(set(A2A6A10A14)-set(A6)-set(A2)-set(A10))
    A11=list(set(A6A7A10A11)-set(A6)-set(A7)-set(A10))
    A1=list(set(A1A2A5A6)-set(A2)-set(A5)-set(A6))
    A3=list(set(A2A3A6A7)-set(A2)-set(A7)-set(A6))

    A5A6A7A8A9A10A11A12=list(set(series8_1_edges))
    A1A2A3A4A5A6A7A8=list(set(series8_2_edges))
    A1A2A5A6A9A10A13A14=list(set(series5_edges))
    A2A3A6A7A10A11A14A15=list(set(series4_edges))

    A4=list(set(A1A2A3A4A5A6A7A8)-set(A1)-set(A2)-set(A3)-set(A5)-set(A6)-set(A7)-set(A8))
    A12=list(set(A5A6A7A8A9A10A11A12)-set(A5)-set(A6)-set(A7)-set(A8)-set(A9)-set(A10)-set(A11))
    A13=list(set(A1A2A5A6A9A10A13A14)-set(A1)-set(A2)-set(A5)-set(A6)-set(A9)-set(A10)-set(A14))
    A15=list(set(A2A3A6A7A10A11A14A15)-set(A2)-set(A3)-set(A6)-set(A7)-set(A10)-set(A11)-set(A14))


    len(set(A13)&set(biopsy_edges))


    ## Export
    with open('../ECv_result/vitroRNAseq/network_th005/viruses_A1.graph.tsv', 'w') as f:
        for i in A1:
            ele=i.split('_')
            row=ele[0]+'\t'+ele[1]+'\n'
            f.write(row)



if __name__ == '__main__':
    main()
