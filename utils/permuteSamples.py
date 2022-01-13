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
meta.to_csv('m.csv', sep=',', index=False)

sample_list = list(meta['Sample'])
#sample_list.sort()

#expr = expr.reindex(columns=['gene'] + sample_list)
expr = expr[['gene'] + sample_list]
expr.to_csv('e.csv', sep=',', index=False)

#sample_id_list_shuffle = np.random.permutation(sample_id_list)
# sort
# sort column names alphabetically
# sort row names alphabetically

# random permute





#expr.columns.get_loc('GLDS_48_Mmus_C57_6J_LVR_GC_C_Rep4_M39_29')
#709
#expr.columns[709]
#'GLDS_48_Mmus_C57_6J_LVR_GC_C_Rep4_M39_29'


