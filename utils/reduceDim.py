import pandas as pd
import sys
import numpy as np
import operator
import argparse


def findMostVaried(df, n):
	# df is genes X samples
	# calculate var, sort cols into n highest vars, drop shape[1]-n cols
	# first find range of var and print to stdout
	if n == 0:
		return df, None
	sdList = df.std(axis=1)
	sdDict = {k: v for v, k in enumerate(sdList)}
	sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=True) 
	topN = sdDictSorted[0:n]
	indices = [x[1] for x in topN]
	slicedDF = df.iloc[indices]
	return slicedDF, indices

def findSumGTDelta(df, delta):
	# first find min sum and print that to stdout
	cSums = df.sum(axis=1)
	cList = list()
	for index, s in cSums.iteritems():
		if s > delta:
			cList.append(index)
	return df.iloc[cList], cList

def removeAlphaZeros(df, alpha):
	# first find max num 0s row and print to stdout
	row_cut_off = int(alpha * len(df.columns))
	df = df[(df == 0).sum(axis='columns') <= row_cut_off]
	return df

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr_file', help='expression file', default=None)
    parser.add_argument('-n', '--num', help='number to reduce to', default=None)
    parser.add_argument('-d', '--delta', help='delta diff b/w highest and lowest expr value', default=None)
    parser.add_argument('-a', '--alpha', help='alpha percentage of 0 expr value', default=None)

    return parser.parse_args()

def main():
	args = parse_args()
	exprFile = args.expr_file
	n = int(args.num)
	delta = int(args.delta)
	alpha = float(int(args.alpha)/100)
	sep=','

	df = pd.read_csv(exprFile, sep=sep, header=0)
	print('original size: ', str(len(df)))

	df_subset, cList = findSumGTDelta(df, delta)
	df_subset = df_subset.reset_index()
	print('after reducing by sum to delta: ', str(delta), str(len(df_subset)))

	df_subset = removeAlphaZeros(df_subset, alpha)
	print('after reducing by removing when percentage zero is at least alpha: ', str(alpha), str(len(df_subset)))

	df_subset,indices = findMostVaried(df_subset, n)
	print('after reducing by n most varied: ', str(n), str(len(df_subset)))

	df_subset = df_subset.drop(columns=['index'])
	genes = df_subset['gene']
	df_subset = df_subset.drop(columns=['gene'])

	df_subset = np.clip(df_subset, 0, a_max=None)

	df_subset.insert(0, 'gene', genes)

	outputFileName = exprFile.split('.csv')[0] + '__reduced_' + str(n) + '_' + str(delta) + '_' + str(sys.argv[4]) + '.csv'

	df_subset.to_csv(outputFileName, sep=',', index=None)

    
if __name__ == "__main__":
	main()	
