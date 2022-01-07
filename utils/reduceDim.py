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

def findSumNear0(df, delta):
	cSums = df.sum(axis=1)
	cList = list()
	for index, s in cSums.iteritems():
		if s > delta:
			cList.append(index)
	return df.iloc[cList], cList
	
def main():
	exprFile = sys.argv[1]
	sep=sys.argv[2]
	n = int(sys.argv[3])
	delta = int(sys.argv[4])

	df = pd.read_csv(exprFile, sep=sep, header=0)

	df_subset, cList = findSumNear0(df, delta)

	df_subset = df_subset.reset_index()

	df_subset,indices = findMostVaried(df_subset, n) 

	df_subset = df_subset.drop(columns=['index'])

	df_subset.to_csv('expr_reduced_' + str(n) + '_' + str(delta) + '.csv', sep=',', index=None)

    
if __name__ == "__main__":
	main()	
