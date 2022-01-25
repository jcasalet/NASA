import pandas as pd
import numpy as np
import argparse
from random import randint

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file anme', default=None, required=True)
    parser.add_argument('-m', '--meta', help='meta file name', default=None, required=True)
    parser.add_argument('-n', '--num', help='number of times more data', default=1)
    parser.add_argument('-v', '--var', help='variance of noise', default=0.1)
    parser.add_argument('-cn', '--colName', help='column name to replicate', default=None)
    parser.add_argument('-cv', '--colVal', help='column value to replicate', default=None)
    return parser.parse_args()

def main():
    options = parse_args()
    expr_df_file = options.expr
    meta_df_file = options.meta
    n = int(options.num)
    var = float(options.var)
    colName = options.colName
    colVal = options.colVal


    expr_df = pd.read_csv(expr_df_file, header=0, sep=',')
    meta_df = pd.read_csv(meta_df_file, header=0, sep=',')
    genes = expr_df['gene']


    if not colName is None:
        col_meta_df = meta_df[meta_df[colName] == colVal]
        # get samples in subset
        col_samples_list = col_meta_df['Sample']
    else:
        list_temp = list(expr_df.columns)
        list_temp.remove('gene')
        col_samples_list = pd.Index(list_temp)
    # subset the expr_df with samples
    expr_samples_df = expr_df[expr_df.columns.intersection(col_samples_list)].T

    expr_df_T = expr_df.T

    # amplify set of samples that match column value
    for i in range(n):
        for sample in col_samples_list:
            # add new sample to expr data
            expr_row = expr_samples_df[expr_samples_df.index == sample]
            noise = np.random.normal(0, var, expr_row.shape)
            noised_expr_row = expr_row + noise
            new_sample = sample + '_' + str(randint(0, 1000000))
            noised_expr_row.rename(index={sample:new_sample}, inplace=True)
            expr_df_T = expr_df_T.append(noised_expr_row, ignore_index=False)

            # add new sample to meta data
            meta_row = meta_df[meta_df['Sample'] == sample]
            new_meta_row = meta_row.copy(deep=True)
            new_meta_row['Sample'] = new_sample
            meta_df = meta_df.append(new_meta_row, ignore_index=True)

    expr_df = expr_df_T.T
    genes = expr_df['gene']
    expr_df = expr_df.drop(columns=['gene'])
    expr_df = np.clip(expr_df, 0, a_max=None)
    expr_df.insert(0, 'gene', genes)
    outputFileNameSuffix = '__expanded_' + str(n) + '_' + str(var) + '.csv'
    expr_df.to_csv(expr_df_file.split('.')[0] + outputFileNameSuffix, index=False)
    meta_df.to_csv(meta_df_file.split('.')[0] + outputFileNameSuffix, index=False)
    print('new expr dims: ', str(expr_df.shape))
    print('new meta dims: ', str(meta_df.shape))



if __name__ == "__main__":
    main()

    '''genes = expr_df['gene']
    expr_df_T = expr_df.T
    expr_df_T_np = expr_df_T.to_numpy()
    df_np = expr_df_T_np[1:]
    orig_df_np = expr_df_T_np[1:]
    orig_info_df = meta_df.copy(deep=True)
    for i in range(0, n):
        noise = np.random.normal(0, var, orig_df_np.shape)
        noised_np = orig_df_np + noise
        noised_np[noised_np<0] = 0
        new_samples = ['sample_' + str(i) + '_' + str(j) for j in range(len(expr_df.columns)-1)]
        orig_info_df['Sample'] = new_samples
        info_df = meta_df.append(orig_info_df)
        df_np = np.concatenate([df_np, noised_np])

    expanded_expr_df = pd.DataFrame(data=df_np, index=meta_df['Sample'], columns=genes).T

    outfilePrefix = 'expanded_' + str(n) + '_' + str(var) + '_'
    expanded_expr_df.to_csv(outfilePrefix + 'expr.csv')
    meta_df.to_csv(outfilePrefix + 'meta.csv', index=None)'''

