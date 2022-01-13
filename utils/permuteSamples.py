import numpy as np
import pandas as pd
import sys

exprFile = sys.argv[1]
metaFile = sys.argv[2]

expr = pd.read_csv(exprFile, sep=',', header=0)
meta = pd.read_csv(metaFile, sep=',', header=0)

sample_id_list = list(range(meta.shape[0]))

sample_id_list_shuffle = np.random.permutation(sample_id_list)

# sort column names alphabetically

# sort row names alphabetically




#expr.columns.get_loc('GLDS_48_Mmus_C57_6J_LVR_GC_C_Rep4_M39_29')
#709
#expr.columns[709]
#'GLDS_48_Mmus_C57_6J_LVR_GC_C_Rep4_M39_29'


