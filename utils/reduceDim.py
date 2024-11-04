import pandas as pd
import sys
import operator
import argparse
from pybiomart import Server
import numpy as np
import matplotlib.pyplot as plt


def transpose_df(df, cur_index_col, new_index_col):
	if not cur_index_col in list(df.columns) or len(df) == 0:
		return df
	df = df.set_index(cur_index_col).T
	df.reset_index(level=0, inplace=True)
	cols = [new_index_col] + list(df.columns)[1:]
	df.columns = cols
	return df

def save_expr(inputDF=None, fileName=None, transpose=False, dropCols=[], cur_index_col=None, new_index_col=None):

	df = inputDF
	if transpose:
		df = transpose_df(df.drop(columns=dropCols), cur_index_col, new_index_col)
		if fileName.endswith('.csv'):
			df.to_csv(fileName, sep=',', index=None)
		elif fileName.endswith('.pkl'):
			df.to_pickle(fileName)
		else:
			print('unknown filename extension: ', fileName)
			sys.exit(1)
	else:
		df = df.drop(columns=dropCols)
		if fileName.endswith('.csv'):
			df.to_csv(fileName, sep=',', index=None)
		elif fileName.endswith('.pkl'):
			df.to_pickle(fileName)
		else:
			print('unknown filename extension: ', fileName)
			sys.exit(1)


def read_expr(fileName, feature_key='gene', sample_key='sample'):
	df = pd.read_csv(fileName, header=0, sep=',')
	first_col = df.columns[0]
	if first_col == feature_key:
		pass
	else:
		df.rename(columns={first_col: feature_key}, inplace=True)

	return transpose_df(df, cur_index_col=feature_key, new_index_col=sample_key)


def callR(cmd):
	import subprocess
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	o, e = proc.communicate(timeout=900)

	print('Output: ' + o.decode('ascii'))
	print('Error: ' + str(e.decode('utf-8')))
	print('code: ' + str(proc.returncode))
	if str(proc.returncode) != '0':
		print('error in callR: exiting')
		sys.exit(1)

def convertIdsToNames(df, outputDir, feature_key='gene', sample_key='sample'):
	input_to_R = outputDir + '/expr_input.csv'
	output_from_R = outputDir + '/expr_output.csv'
	save_expr(inputDF=df, fileName=input_to_R, transpose=True, dropCols=[], cur_index_col=sample_key,
				   new_index_col=feature_key)
	R_cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/convert_id_to_gene.R', '--args', input_to_R, output_from_R]
	callR(R_cmd)
	return read_expr(output_from_R)


def filterGenesByType(df, sample_key, feature_key, gene_type='protein_coding', id='id'):
	gene_types = {'ribozyme', 'protein_coding', 'rRNA', 'TEC', 'IG_D_pseudogene', 'snRNA', 'IG_LV_gene', 'pseudogene',
				  'IG_J_gene', 'transcribed_unitary_pseudogene', 'processed_pseudogene', 'IG_V_gene', 'Mt_tRNA',
				  'TR_J_pseudogene', 'miRNA', 'Mt_rRNA', 'sRNA', 'IG_C_pseudogene', 'IG_C_gene', 'TR_J_gene',
				  'IG_pseudogene', 'transcribed_processed_pseudogene', 'scRNA', 'lncRNA', 'TR_V_pseudogene',
				  'TR_V_gene', 'misc_RNA', 'TR_D_gene', 'translated_unprocessed_pseudogene',
				  'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene',
				  'IG_V_pseudogene', 'scaRNA', 'TR_C_gene', 'IG_D_gene', 'snoRNA'}
	if not gene_type in gene_types:
		print('gene_type: ' + str(gene_type) + ' not recognized')
		sys.exit(1)
	server = Server(host='http://www.ensembl.org')
	dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl'])
	gene_info = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'])

	if id == 'id':
		filter_genes = list(
			gene_info[(gene_info['Gene type'] == gene_type) & (gene_info['Gene type'] != 'Mt_rRNA')]['Gene stable ID'])
		filter_columns = [sample_key] + filter_genes
		df = df[df.columns.intersection(filter_columns)]
		new_columns = list(df.drop(columns=[sample_key]))
		# gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene name'])
		gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene stable ID'])
		df.columns = [sample_key] + gene_names
	elif id == feature_key:
		filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene name'])
		filter_columns = [sample_key] + filter_genes
		df = df[df.columns.intersection(filter_columns)]
	df = df.loc[:, df.columns.notna()]
	return df

def filter_cvs(df, coef_var=0.5):
	# df is genes X samples
	# calculate coefficient of variation
	cvs=list()
	for i in range(len(df)):
		m=np.mean(df.iloc[i][1:])
		sd=np.std(df.iloc[i][1:])
		cvs.append(sd/m)

	'''# plot hist of dist of coev of variation
	fig, axs = plt.subplots()
	axs.hist(cvs, bins=20)
	plt.show()'''

	# keep genes with cv > coef_var
	indices = list()
	for i in range(len(cvs)):
		if cvs[i] > coef_var:
			indices.append(i)

	return df.iloc[indices]

def findMostVaried(df, n, key):
	# df is genes X samples
	# calculate var, sort cols into n highest vars, drop shape[1]-n cols
	# first find range of var and print to stdout
	if n == 0:
		return df, None
	if 'index' in df.columns:
		df.reset_index(drop=True, inplace=True)
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
	df.reset_index(drop=True, inplace=True)
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
	parser.add_argument('-c', '--coef_var', help='coefficient of variation threshold', default=0.0)
	parser.add_argument('-e', '--expr_file', help='expression file', default=None)
	parser.add_argument('-n', '--num', help='number to reduce to', default=0)
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
	coef_var = float(args.coef_var)
	sep=','
	key = args.key

	df = pd.read_csv(exprFile, sep=sep, header=0)
	print('original size: ', str(df.shape))

	df = convertIdsToNames(df=df, outputDir='/tmp')
	print('dims after converting to names: ', df.shape)

	df = filterGenesByType(df, key, feature_key='gene', gene_type='protein_coding', id='gene')
	print('after reducing by removing non-protein-coding genes: ', df.shape)

	df = transpose_df(df, cur_index_col='sample', new_index_col='gene')

	df = removeAlphaZeros(df, alpha)
	print('after reducing by removing when percentage zero is at least alpha: ', str(alpha), df.shape)

	df, cList = findSumGTSigma(df, sigma)
	print('after reducing by sum to sigma: ', str(sigma), df.shape)

	df = removeDeltaDiff(df, delta)
	print('after reducing by removing when (max - min) is at most delta: ', str(delta), df.shape)

	df,indices = findMostVaried(df, n, key)
	print('after reducing by n most varied: ', str(n), df.shape)

	df = filter_cvs(df, coef_var)
	print('after reducing by coef_var : ', str(coef_var), df.shape)

	if 'index' in list(df.columns):
		df.drop(columns=['index'], inplace=True)
	outputFileName = exprFile.split('.csv')[0] + "__reduced_" +  \
					 "_a_" + str(alpha) + \
					 "_s_" + str(sigma) +  \
					 "_d_" + str(delta) + \
					 "_n_" + str(n)  + \
					 ".csv"

	df.to_csv(outputFileName, sep=',', index=None)

    
if __name__ == "__main__":
	main()	
