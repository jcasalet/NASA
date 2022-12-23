from synthetic.my_synthetic_generator import synthetic_generator

N_mc = 1

for n in [5,10,50]:
    for i in range(N_mc):
        synthetic_df, colors = synthetic_generator(n=30,d_layer=3,n_layer=[5,10,20],mu=0,sigma=1,n_causal=n)
        synthetic_name = 'full_fw_synthetic_sem_n_causal_'+str(n)+'_'+str(i)+'.pickle'
        synthetic_loc = 'data/synthetic/'+synthetic_name
        synthetic_df.to_pickle(synthetic_loc)

print("Generated synthetic datasets")

import matplotlib.pyplot as plt
import numpy as np

plotName='bob'
crispGenesFileName='features.txt'
envVar = 'env_split'
sampleName = 'Subj_ID'
targetName = 'Target'
logBase = '2'
includeText = False

df=synthetic_df
means_list_dict=dict()
vars_list_dict=dict()
means_dict=dict()
vars_dict=dict()

for env in set(df[envVar]):
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

	if logBase is None:
		plt.scatter(x=means_list_dict[env], y=vars_list_dict[env])
	elif int(logBase) == 10:
		plt.scatter(x=np.log10([i+1 for i in means_list_dict[env]]), y=np.log10([j+1 for j in vars_list_dict[env]]))
	elif int(logBase) == 2:
		plt.scatter(x=np.log2([i+1 for i in means_list_dict[env]]), y=np.log2([j+1 for j in vars_list_dict[env]]))

	crisp_x = list()
	crisp_y = list()

	for gene in crispGenes:
		if logBase is None:
			crisp_x.append(means_dict[env][gene])
			crisp_y.append(vars_dict[env][gene])
			x = means_dict[env][gene]
			y = vars_dict[env][gene]
		elif int(logBase) == 10:
			crisp_x.append(np.log10(1+means_dict[env][gene]))
			crisp_y.append(np.log10(1+vars_dict[env][gene]))
			x = np.log10(1+means_dict[env][gene])
			y = np.log10(1+vars_dict[env][gene])
		elif int(logBase) == 2:
			crisp_x.append(np.log2(1+means_dict[env][gene]))
			crisp_y.append(np.log2(1+vars_dict[env][gene]))
			x = np.log2(1+means_dict[env][gene])
			y = np.log2(1+vars_dict[env][gene])

		if bool(includeText):
			plt.text(x, y, gene, color='black')

	plt.scatter(x=crisp_x, y=crisp_y, marker='*', color='red', s=100)


	plt.savefig(plotName + '_' + str(env) + '.png')
	plt.close()

nrows = colors.shape[0]
ncols = colors.shape[1]
print('rows = ', nrows)
print('cols = ', ncols)
for i in range(colors.shape[0]):
	for j in range(colors.shape[1]):
		if colors[i][j] == 0:
			c='red'
		elif colors[i][j] == 1:
			c='blue'
		elif colors[i][j] == 2:
			c='orange'
		elif colors[i][j] == 3:
			c='purple'
		elif colors[i][j] == 4:
			c='brown'
		else:
			c='black'
		plt.scatter(x=i, y=j, marker='.', color=c, s=100)
plt.savefig(plotName + '_' + 'colors.png')
plt.close()
#####################

