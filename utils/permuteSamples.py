import numpy as np
import pandas as pd
import sys

exprFile = sys.argv[1]
metaFile = sys.argv[2]

expr = pd.read_csv(exprFile, sep=',', header=0)
meta = pd.read_csv(metaFile, sep=',', header=0)

sample2index_dict = dict()

for i in range(meta.shape[0]):
    sample = meta.iloc[i]['Sample']
    j = expr.columns.get_loc(sample)
    sample2index_dict[sample] = (i, j)



meta = meta.sample(frac=1)
meta.to_csv(metaFile.split('.csv')[0] + '_permuted.csv', sep=',', index=False)

sample_list = list(meta['Sample'])
#sample_list.sort()

#expr = expr.reindex(columns=['gene'] + sample_list)
expr = expr[['gene'] + sample_list]
expr.to_csv(exprFile.split('.csv')[0] + '_permuted.csv', sep=',', index=False)



