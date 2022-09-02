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

if 'env' in df.columns: 
	df.drop(columns=['env'], inplace=True)

if 'oro_thresh' in df.columns: 
	df.drop(columns=['oro_thresh'], inplace=True)

means_dict = dict(df.mean())
vars_dict = dict(df.var())

means_list = list(means_dict.values())
vars_list = list(vars_dict.values())
#####################


#####################
with open(crispGenesFileName, 'r') as f:
	crispGenes = f.read().splitlines()
f.close()
print(crispGenes)
#####################


#####################
fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel('mean', fontsize=18)
ax.set_ylabel('variance', fontsize=18)
fig.suptitle(plotName, fontsize=24)

plt.scatter(x=np.log10([i+1 for i in means_list]), y=np.log10([j+1 for j in vars_list]))
#plt.scatter(x=means_list, y=vars_list)

crisp_x = list()
crisp_y = list()

print(means_dict)

for gene in crispGenes:
	#crisp_x.append(np.log10(1+means_dict[gene]))
	#crisp_x.append(means_dict[gene])
	#crisp_y.append(np.log10(1+vars_dict[gene]))
	#crisp_y.append(vars_dict[gene])
	x = np.log10(1+means_dict[gene])
	y = np.log10(1+vars_dict[gene])
	plt.text(x, y, gene)

#plt.scatter(x=crisp_x, y=crisp_y, marker='*', color='red', s=100)

plt.savefig(plotName + '.png')
#####################
