import numpy as np
import pandas as pd
import sys

expr_df = pd.read_csv(sys.argv[1], sep=',', header=0, index_col=0)
print('expr_df shape: ', str(expr_df.shape))
meta_df = pd.read_csv(sys.argv[2], sep=',', header=0, index_col=0)
print('meta_df shape: ', str(meta_df.shape))
x_gen = np.genfromtxt(sys.argv[3], delimiter=',')
print('x_gen shape: ', str(x_gen.shape))
num_test = int(sys.argv[4])

expr_df_genes = expr_df.index
index = meta_df.shape[0] - num_test
expr_df_samples=meta_df.index[index:]
x_gen_df = pd.DataFrame(data=x_gen.T, index=expr_df_genes, columns=expr_df_samples)

x_gen_df.to_csv(sys.argv[3].split('.csv')[0] + '_pd.csv', sep=',', header=True, index=True)
meta_df[index:].to_csv('m_gen_pd.csv', sep=',', header=True, index=True)