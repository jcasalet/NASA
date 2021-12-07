import pandas as pd
import numpy as np
import sys

expr_df_file = sys.argv[1]
info_df_file = sys.argv[2]
n = int(sys.argv[3])
var = float(sys.argv[4])

expr_df = pd.read_csv(expr_df_file, header=0, sep=',')
info_df = pd.read_csv(info_df_file, header=0, sep=',')
# need to randomize the selection of samples?!

genes = expr_df['gene']
expr_df_T = expr_df.T
expr_df_T_np = expr_df_T.to_numpy()
df_np = expr_df_T_np[1:]
orig_df_np = expr_df_T_np[1:]
orig_info_df = info_df.copy(deep=True)
for i in range(0, n):
    noise = np.random.normal(0, var, orig_df_np.shape)
    noised_np = orig_df_np + noise
    noised_np[noised_np<0] = 0
    new_samples = ['sample_' + str(i) + '_' + str(j) for j in range(len(expr_df.columns)-1)]
    orig_info_df['Sample'] = new_samples
    info_df = info_df.append(orig_info_df)
    df_np = np.concatenate([df_np, noised_np])

expanded_expr_df = pd.DataFrame(data=df_np, index=info_df['Sample'], columns=genes).T

expanded_expr_df.to_csv('expanded_expr_df.csv')
info_df.to_csv('expanded_info_df.csv', index=None)

