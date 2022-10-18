import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

dfFileName=sys.argv[1]
plotName=sys.argv[2]
crispGenesFileName=sys.argv[3]


#####################
with open(dfFileName, 'rb') as f:
	#df=pd.read_csv(f, header=0, sep=',')
	df=pd.read_pickle(f)
f.close()
means_list_dict=dict()
vars_list_dict=dict()
means_dict=dict()
vars_dict=dict()
for env in df['env']:
	df_env = df[df['env'] == env]
	if 'env' in df_env.columns:
		df_env.drop(columns=['env'], inplace=True)

	if 'oro_thresh' in df_env.columns:
		df_env.drop(columns=['oro_thresh'], inplace=True)

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
for env in df['env']:
	fig, ax = plt.subplots(figsize=(10,10))
	ax.set_xlabel('mean', fontsize=18)
	ax.set_ylabel('variance', fontsize=18)
	fig.suptitle(plotName, fontsize=24)

	plt.scatter(x=np.log10([i+1 for i in means_list_dict[env]]), y=np.log10([j+1 for j in vars_list_dict[env]]))
	#plt.scatter(x=means_list_dict[env], y=vars_list_dict[env])
	#plt.scatter(x=[i+1 for i in means_list_dict[env]], y=[j+1 for j in vars_list_dict[env]])


	crisp_x = list()
	crisp_y = list()

	for gene in crispGenes:
		crisp_x.append(np.log10(1+means_dict[env][gene]))
		#crisp_x.append(means_dict[env][gene])
		crisp_y.append(np.log10(1+vars_dict[env][gene]))
		#crisp_y.append(vars_dict[env][gene])
		#x = np.log10(1+means_dict[env][gene])
		#y = np.log10(1+vars_dict[env][gene])
		#plt.text(x, y, gene, color='red', )

	#plt.scatter(x=crisp_x, y=crisp_y, marker='*', color='red', s=100)


	plt.savefig(plotName + '_' + env + '.png')
#####################
