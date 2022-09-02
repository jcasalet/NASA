import pandas as pd
import sys
import scipy.stats as stats
import statistics
import matplotlib.pyplot as plt
import pylab

sigThresh = 0.05

#expr=pd.read_csv('liver_env=mission:dissection:library-prep.csv', sep=',', header=0)
expr=pd.read_pickle(sys.argv[1])
meta=pd.read_csv(sys.argv[2], sep=',', header=0)
with open(sys.argv[3], 'r') as f:
	genes=f.read().splitlines()
f.close()

flight_samples=list(meta[meta['group']=='Flight']['sample'])
ground_samples=list(meta[meta['group']=='Ground']['sample'])
basal_samples=list(meta[meta['group']=='Basal']['sample'])
vivarium_samples=list(meta[meta['group']=='Vivarium']['sample'])

expr_per_gene_per_group=dict()
results = dict()

for gene in genes:
	expr_per_gene_per_group[gene] = dict()
	expr_per_gene_per_group[gene]['flight'] = list()
	expr_per_gene_per_group[gene]['ground'] = list()
	expr_per_gene_per_group[gene]['basal'] = list()
	expr_per_gene_per_group[gene]['vivarium'] = list()
	for sample in expr['sample']:
		num = float(expr[expr['sample']==sample][gene])
		if sample in flight_samples:
			expr_per_gene_per_group[gene]['flight'].append(num)
		elif sample in ground_samples:
			expr_per_gene_per_group[gene]['ground'].append(num)
		elif sample in basal_samples:
			expr_per_gene_per_group[gene]['basal'].append(num)
		elif sample in vivarium_samples:
			expr_per_gene_per_group[gene]['vivarium'].append(num)
		else:
			print('unknown sample: ', sample)

	results[gene] = dict()
	for group in ['flight', 'ground', 'basal', 'vivarium']:
		results[gene][group] = dict()
		results[gene][group]['mean'] = float('%.3f'%(statistics.mean(expr_per_gene_per_group[gene][group])))
		results[gene][group]['variance'] = float('%.3f'%(statistics.variance(expr_per_gene_per_group[gene][group])))
		results[gene][group]['n'] = len(expr_per_gene_per_group[gene][group])


	nonflight = expr_per_gene_per_group[gene]['ground'] + expr_per_gene_per_group[gene]['basal'] + expr_per_gene_per_group[gene]['vivarium']
	results[gene]['nonflight'] = dict()
	results[gene]['nonflight']['mean'] = float('%.3f'%(statistics.mean(nonflight)))
	results[gene]['nonflight']['variance'] = float('%.3f'%(statistics.variance(nonflight)))
	results[gene]['nonflight']['n'] = len(nonflight)

	'''flight_ground_t_test = (flight_mean - ground_mean) / (math.sqrt(flight_var/n_flight + ground_var/n_ground))
	flight_basal_t_test = (flight_mean - basal_mean) / (math.sqrt(flight_var/n_flight + basal_var/n_basal))
	flight_vivarium_t_test = (flight_mean - vivarium_mean) / (math.sqrt(flight_var/n_flight + vivarium_var/n_vivarium))'''
	results[gene]['flight_vs_ground'] = dict()
	results[gene]['flight_vs_ground']['t-test'] = float('%.3f'%(stats.ttest_ind(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['ground']).pvalue))
	results[gene]['flight_vs_ground']['wilcoxon'] = float('%.3f'%(stats.ranksums(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['ground']).pvalue))

	results[gene]['flight_vs_basal'] = dict()
	results[gene]['flight_vs_basal']['t-test'] = float('%.3f'%(stats.ttest_ind(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['basal']).pvalue))
	results[gene]['flight_vs_basal']['wilcoxon'] = float('%.3f'%(stats.ranksums(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['basal']).pvalue))

	results[gene]['flight_vs_vivarium'] = dict()
	results[gene]['flight_vs_vivarium']['t-test'] = float('%.3f'%(stats.ttest_ind(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['vivarium']).pvalue))
	results[gene]['flight_vs_vivarium']['wilcoxon'] = float('%.3f'%(stats.ranksums(expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['vivarium']).pvalue))

	results[gene]['flight_vs_nonflight'] = dict()
	results[gene]['flight_vs_nonflight']['t-test'] = float('%.3f'%(stats.ttest_ind(expr_per_gene_per_group[gene]['flight'], nonflight).pvalue))
	results[gene]['flight_vs_nonflight']['wilcoxon'] = float('%.3f'%(stats.ranksums(expr_per_gene_per_group[gene]['flight'], nonflight).pvalue))

	'''print('gene: ', gene)

	print('flight - ground = ', results[gene]['flight']['mean'] - results[gene]['ground']['mean'])
	print('flight/ground welch t-test: ',  results[gene]['flight_vs_ground']['t-test'])
	print('flight/ground wilcoxon rank sum t-test: ', results[gene]['flight_vs_ground']['wilcoxon'])


	print('flight - basal = ', results[gene]['flight']['mean'] - results[gene]['basal']['mean'])
	print('flight/basal welch t-test: ',  results[gene]['flight_vs_basal']['t-test'])
	print('flight/basal wilcoxon rank sum t-test: ',  results[gene]['flight_vs_basal']['wilcoxon'])


	print('flight - vivarium = ', results[gene]['flight']['mean'] - results[gene]['vivarium']['mean'])
	print('flight/vivarium welch t-test: ',  results[gene]['flight_vs_vivarium']['t-test'])
	print('flight/vivarium wilcoxon rank sum t-test: ',  results[gene]['flight_vs_vivarium']['wilcoxon'])

	print('flight - nonflight = ', results[gene]['flight']['mean'] - results[gene]['nonflight']['mean'])
	print('flight/nonflight welch t-test: ',  results[gene]['flight_vs_nonflight']['t-test'])
	print('flight/nonflight wilcoxon rank sum t-test: ',  results[gene]['flight_vs_nonflight']['wilcoxon'])'''



	nbins = 10
	'''fig, axs = plt.subplots(1, 3)
	axs[0].hist(expr_per_gene_per_group[gene]['flight'], bins=nbins)
	axs[0].set_title('flight')
	axs[1].hist(expr_per_gene_per_group[gene]['ground'], bins=nbins)
	axs[1].set_title('ground')
	axs[2].hist(nonflight, bins=nbins)
	axs[2].set_title('non-flight')
	fig.suptitle(gene)
	plt.savefig(gene + '.png', dpi=300)'''

	'''colors = ['red', 'blue', 'orange']
	combo_data = [expr_per_gene_per_group[gene]['flight'], expr_per_gene_per_group[gene]['ground'], nonflight]
	plt.hist(x=combo_data, bins=nbins, color=colors, density=True, stacked=True)
	plt.legend(['flight', 'ground', 'non-flight'])
	plt.title(gene)
	plt.savefig(gene + '.png', dpi=300)'''


	# qqplot
	'''stats.probplot(expr_per_gene_per_group[gene]['flight'])
	stats.probplot(expr_per_gene_per_group[gene]['ground'])
	stats.probplot(nonflight)
	pylab.show()'''

import json
print(json.dumps(results, sort_keys=False, indent=4))