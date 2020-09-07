
import pandas as pd

def main():
    ecv5=pd.read_table('deltaECv_table_Series5_th1.txt', sep='\t', header=0)

    with open('test.raph.tsv', 'w') as f:
        for i in ecv1_edges:
            row = i[0] + '\t' + i[1] + '\n'
            f.write(row)



if __name__ = '__main__':
    main()

