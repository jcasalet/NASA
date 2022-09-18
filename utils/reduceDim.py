import pandas as pd
import sys
import numpy as np
import operator
import argparse


def findMostVaried(df, n, key):
	# df is genes X samples
	# calculate var, sort cols into n highest vars, drop shape[1]-n cols
	# first find range of var and print to stdout
	if n == 0:
		return df, None
	index = df.index
	sdList = df.std(axis=1)
	sdDict = {k: v for v, k in enumerate(sdList)}
	sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=True) 
	topN = sdDictSorted[0:n]
	indices = [x[1] for x in topN]
	slicedDF = df.iloc[indices]
	return slicedDF, indices

def findSumGTSigma(df, sigma):
	if sigma == 0:
		return df, None
	# first find min sum and print that to stdout
	df.reset_index(inplace=True)
	cSums = df.sum(axis=1)
	cList = list()
	for index, s in cSums.iteritems():
		if s > sigma:
			cList.append(index)
	temp = df.iloc[cList]
	return temp, cList

def removeAlphaZeros(df, alpha):
	if alpha == 0:
		return df
	return df[(df == 0).sum(axis='columns') <= int(alpha * len(df.columns))]

def removeDeltaDiff(df, delta):
	if delta == 0:
		return df
	return df[df.max(axis=1) - df.min(axis=1) > delta]

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-e', '--expr_file', help='expression file', default=None)
	parser.add_argument('-n', '--num', help='number to reduce to', default=None)
	parser.add_argument('-d', '--delta', help='delta diff of expr vals max to min across samples', default=0)
	parser.add_argument('-s', '--sigma', help='sigma sum of expr vals across samples', default=0)
	parser.add_argument('-a', '--alpha', help='alpha percentage of 0 expr value', default=0)
	parser.add_argument('-k', '--key', help='name of key column', default='gene')
	return parser.parse_args()

def main():
	args = parse_args()
	exprFile = args.expr_file
	n = int(args.num)
	delta = int(args.delta)
	alpha = float(int(args.alpha)/100)
	sigma = int(args.sigma)
	sep=','
	key = args.key

	df = pd.read_csv(exprFile, sep=sep, header=0)
	print('original size: ', str(df.shape))

	df = removeAlphaZeros(df, alpha)
	print('after reducing by removing when percentage zero is at least alpha: ', str(alpha), str(len(df)))

	df, cList = findSumGTSigma(df, sigma)
	print('after reducing by sum to sigma: ', str(sigma), str(len(df)))

	df = removeDeltaDiff(df, delta)
	print('after reducing by removing when (max - min) is at most delta: ', str(delta), str(len(df)))

	df,indices = findMostVaried(df, n, key)
	print('after reducing by n most varied: ', str(n), str(len(df)))

	if 'index' in list(df.columns):
		df.drop(columns=['index'], inplace=True)
	outputFileName = exprFile.split('.csv')[0] + '__reduced_' + str(n) + '_' + str(delta) + '_' + str(alpha) + '.csv'
	df.to_csv(outputFileName, sep=',', index=None)

    
if __name__ == "__main__":
	main()	
