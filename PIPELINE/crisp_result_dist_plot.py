import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

#envVar = 'env_split'
#sampleName = 'Subj_ID'
#targetName = 'Target'

parser = argparse.ArgumentParser()
parser.add_argument('-e', '--env', help='env field', default=None, required=True)
parser.add_argument('-s', '--sample', help='sample field', default=None, required=True)
parser.add_argument('-t', '--target', help='target field', default=None, required=True)
parser.add_argument('-df', '--dataFile', help='data file name', default=None, required=True)
parser.add_argument('-gf', '--geneFile', help='gene file name', default=None, required=True)
parser.add_argument('-pf', '--plotFile', help='plot file name', default=None, required=True)
parser.add_argument('-l', '--logBase', help='log base', default=None, required=False)
parser.add_argument('-i', '--includeText', help='include text or not', default=False, required=False)


args = parser.parse_args()

dfFileName=args.dataFile
plotName=args.plotFile
crispGenesFileName=args.geneFile
envVar = args.env
sampleName = args.sample
targetName = args.target



#####################
with open(dfFileName, 'rb') as f:
	#df=pd.read_csv(f, header=0, sep=',')
	df=pd.read_pickle(f)
f.close()
means_list_dict=dict()
vars_list_dict=dict()
means_dict=dict()
vars_dict=dict()
for env in df[envVar]:
	df_env = df[df[envVar] == env]
	if envVar in df_env.columns:
		df_env.drop(columns=[envVar], inplace=True)

	if targetName in df_env.columns:
		df_env.drop(columns=[targetName], inplace=True)

	if sampleName in df_env.columns:
		df_env.drop(columns=[sampleName], inplace=True)

	means_dict[env] = dict(df_env.mean())
	vars_dict[env] = dict(df_env.var())

	means_list_dict[env] = list(means_dict[env].values())
	vars_list_dict[env] = list(vars_dict[env].values())



#####################


#####################
with open(crispGenesFileName, 'r') as f:
	crispGenes = f.read().splitlines()
f.close()
print(crispGenes)
#####################




#####################
for env in df[envVar]:
	fig, ax = plt.subplots(figsize=(10,10))
	ax.set_xlabel('mean', fontsize=18)
	ax.set_ylabel('variance', fontsize=18)
	fig.suptitle(plotName, fontsize=24)

	if args.logBase is None:
		plt.scatter(x=means_list_dict[env], y=vars_list_dict[env])
	elif int(args.logBase) == 10:
		plt.scatter(x=np.log10([i+1 for i in means_list_dict[env]]), y=np.log10([j+1 for j in vars_list_dict[env]]))
	elif int(args.logBase) == 2:
		plt.scatter(x=np.log2([i+1 for i in means_list_dict[env]]), y=np.log2([j+1 for j in vars_list_dict[env]]))

	crisp_x = list()
	crisp_y = list()

	for gene in crispGenes:
		if args.logBase is None:
			crisp_x.append(means_dict[env][gene])
			crisp_y.append(vars_dict[env][gene])
			x = means_dict[env][gene]
			y = vars_dict[env][gene]
		elif int(args.logBase) == 10:
			crisp_x.append(np.log10(1+means_dict[env][gene]))
			crisp_y.append(np.log10(1+vars_dict[env][gene]))
			x = np.log10(1+means_dict[env][gene])
			y = np.log10(1+vars_dict[env][gene])
		elif int(args.logBase) == 2:
			crisp_x.append(np.log2(1+means_dict[env][gene]))
			crisp_y.append(np.log2(1+vars_dict[env][gene]))
			x = np.log2(1+means_dict[env][gene])
			y = np.log2(1+vars_dict[env][gene])

		if bool(args.includeText):
			plt.text(x, y, gene, color='black')

	plt.scatter(x=crisp_x, y=crisp_y, marker='*', color='red', s=100)


	plt.savefig(plotName + '_' + str(env) + '.png')
#####################
