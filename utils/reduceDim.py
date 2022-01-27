import pandas as pd
import sys
import numpy as np
import operator


def findMostVaried(df, n):
	# df is genes X samples
	# calculate var, sort cols into n highest vars, drop shape[1]-n cols
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
	cSums = df.sum(axis=1)
	cList = list()
	for index, s in cSums.iteritems():
		if s > delta:
			cList.append(index)
	return df.iloc[cList], cList

def removeAlphaZeros(df, alpha):
	row_cut_off = int(alpha * len(df.columns))
	df = df[(df == 0).sum(axis='columns') <= row_cut_off]
	return df


def main():
	exprFile = sys.argv[1]
	n = int(sys.argv[2])
	delta = int(sys.argv[3])
	alpha = float(sys.argv[4])
	sep=','

	df = pd.read_csv(exprFile, sep=sep, header=0)

	df_subset, cList = findSumGTDelta(df, delta)

	df_subset = df_subset.reset_index()

	print(len(df_subset))
	df_subset = removeAlphaZeros(df_subset, alpha)
	print(len(df_subset))

	df_subset,indices = findMostVaried(df_subset, n)

	df_subset = df_subset.drop(columns=['index'])

	genes = df_subset['gene']

	df_subset = df_subset.drop(columns=['gene'])

	df_subset = np.clip(df_subset, 0, a_max=None)

	df_subset.insert(0, 'gene', genes)

	outputFileName = exprFile.split('.')[0] + '__reduced_' + str(n) + '_' + str(delta) + '.csv'

	df_subset.to_csv(outputFileName, sep=',', index=None)

    
if __name__ == "__main__":
	main()	
