import pandas as pd
import numpy as np
import argparse
from random import randint
from multiprocessing import Process, Queue, cpu_count


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file anme', default=None, required=True)
    parser.add_argument('-m', '--meta', help='meta file name', default=None, required=True)
    parser.add_argument('-n', '--num', help='number of times more data', default=1)
    parser.add_argument('-v', '--var', help='variance of noise', default=0.1)
    parser.add_argument('-cn', '--colName', help='column name to replicate', default=None)
    parser.add_argument('-cv', '--colVal', help='column value to replicate', default=None)
    parser.add_argument("-nc", "--ncpu", help="Number of processes. Default=cpu_count()", default=cpu_count())

    return parser.parse_args()

def divide(n, d):
   res = list()
   qu = int(n/d)
   rm = n%d
   for i in range(d):
       if i < rm:
           res.append(qu + 1)
       else:
           res.append(qu)
   return res

def getStartAndEnd(partitionSizes, threadID):
    start = 0
    for i in range(threadID):
        start += partitionSizes[i]
    end = start + partitionSizes[threadID]

    return start, end

def process_samples_subset(q, samples, expr_samples_df, expr_df_T, meta_df, n, var, threadID, numProcs):

    results = dict()
    partitionSizes = divide(len(samples), numProcs)
    start, end = getStartAndEnd(partitionSizes, threadID)
    print('id: ', str(threadID), 'numProcs: ', str(numProcs), 'start: ', str(start), 'end: ', str(end))
    temp_meta_df = pd.DataFrame(columns=meta_df.columns)
    temp_expr_df = pd.DataFrame(columns=expr_df_T.columns)
    for i in range(n):
        for sample in samples[start:end]:
            # add new sample to expr data
            expr_row = expr_samples_df[expr_samples_df.index == sample]
            noise = np.random.normal(0, var, expr_row.shape)
            noised_expr_row = expr_row + noise
            new_sample = sample + '_' + str(randint(0, 1000000))
            #new_sample = sample + '_' + str((i+1) * n)
            noised_expr_row.rename(index={sample: new_sample}, inplace=True)
            temp_expr_df = temp_expr_df.append(noised_expr_row, ignore_index=False)

            # add new sample to meta data
            meta_row = meta_df[meta_df['Sample'] == sample]
            new_meta_row = meta_row.copy(deep=True)
            new_meta_row['Sample'] = new_sample
            temp_meta_df = temp_meta_df.append(new_meta_row, ignore_index=True)
    #return expr_df_T, meta_df
    results['expr_df_T'] = temp_expr_df
    results['meta_df'] = temp_meta_df
    print('subset expr dims: ', str(temp_expr_df.shape))
    print('subset meta dims: ', str(temp_meta_df.shape))
    q.put(results)

def main():
    options = parse_args()
    expr_df_file = options.expr
    meta_df_file = options.meta
    n = int(options.num)
    var = float(options.var)
    colName = options.colName
    colVal = options.colVal

    numProcs = int(options.ncpu)

    expr_df = pd.read_csv(expr_df_file, header=0, sep=',')
    meta_df = pd.read_csv(meta_df_file, header=0, sep=',')

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

    # save gene list
    genes = expr_df['gene']
    # transpose expr_df to have rows be samples
    expr_df_T = expr_df.T

    # amplify set of samples that match column value
    q = Queue()
    processList = list()
    for i in range(numProcs):
        p = Process(target=process_samples_subset, args=(q,col_samples_list, expr_samples_df, expr_df_T, meta_df, n,
                                                         var, i, numProcs, ))
        #expr_df_T, meta_df = process_samples_subset(col_samples_list, expr_samples_df, expr_df_T, meta_df, var)
        p.start()
        processList.append(p)

    results = dict()
    for i in range(numProcs):
        results.update(q.get())
        print('results expr shape: ', str(results['expr_df_T'].shape))
        print('results meta shape: ', str(results['meta_df'].shape))

        expr_df_T = expr_df_T.append(results['expr_df_T'], ignore_index=False )
        meta_df = meta_df.append(results['meta_df'], ignore_index=False)

        print('expr_df_T after append dims:', str(expr_df_T.shape))
        print('meta_df after append dims: ', str(meta_df.shape))

    for i in range(numProcs):
        print('joing thread: ', str(i))
        processList[i].join()

    expr_df = expr_df_T.T
    genes = expr_df['gene']
    expr_df = expr_df.drop(columns=['gene'])
    expr_df = np.clip(expr_df, 0, a_max=None)
    expr_df.insert(0, 'gene', genes)
    outputFileNameSuffix = '__expanded_' + str(n) + '_' + str(var) + '_' + str(colVal) + '.csv'
    expr_df.to_csv(expr_df_file.split('.')[0] + outputFileNameSuffix, index=False)
    meta_df.to_csv(meta_df_file.split('.')[0] + outputFileNameSuffix, index=False)
    print('new expr dims: ', str(expr_df.shape))
    print('new meta dims: ', str(meta_df.shape))


if __name__ == "__main__":
    main()
