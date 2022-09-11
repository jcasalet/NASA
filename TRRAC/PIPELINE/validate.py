import pandas as pd
import sys
import scipy.stats as stats
import statistics
import matplotlib.pyplot as plt
import pylab
from sklearn.preprocessing import MinMaxScaler
import numpy as np

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
scaler = MinMaxScaler(feature_range=(0,1))

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
		counts_2d = []
		for c in expr_per_gene_per_group[gene][group]:
			counts_2d.append([c])
		results[gene][group]['normalized'] = list(scaler.fit_transform(counts_2d).flatten())
		#results[gene][group]['mean'] = float('%.3f'%(statistics.mean(expr_per_gene_per_group[gene][group])))
		#results[gene][group]['variance'] = float('%.3f'%(statistics.variance(expr_per_gene_per_group[gene][group])))
		results[gene][group]['mean'] = float('%.3f' % (statistics.mean(results[gene][group]['normalized'])))
		results[gene][group]['variance'] = float('%.3f' % (statistics.variance(results[gene][group]['normalized'])))
		results[gene][group]['n'] = len(expr_per_gene_per_group[gene][group])


	nonflight = expr_per_gene_per_group[gene]['ground'] + expr_per_gene_per_group[gene]['basal'] + expr_per_gene_per_group[gene]['vivarium']
	counts_2d = []
	for c in nonflight:
		counts_2d.append([c])
	results[gene]['nonflight'] = dict()
	results[gene]['nonflight']['normalized'] = list(scaler.fit_transform(counts_2d).flatten())
	results[gene]['nonflight']['mean'] = float('%.3f'%(statistics.mean(results[gene]['nonflight']['normalized'])))
	results[gene]['nonflight']['variance'] = float('%.3f'%(statistics.variance(results[gene]['nonflight']['normalized'])))
	results[gene]['nonflight']['n'] = len(nonflight)


	for group in ['ground', 'basal', 'vivarium', 'nonflight']:
		comparison = 'flight_vs_' + group
		results[gene][comparison] = dict()
		results[gene][comparison]['t-test'] = float('%.3f' % (stats.ttest_ind(results[gene]['flight']['normalized'], results[gene][group]['normalized']).pvalue))
		results[gene][comparison]['wilcoxon'] = float('%.3f' % (stats.ranksums(results[gene]['flight']['normalized'], results[gene][group]['normalized']).pvalue))
		results[gene][comparison]['ks-test'] = float('%.3f' % (stats.kstest(np.log(results[gene]['flight']['normalized']), np.log(results[gene][group]['normalized'])).pvalue))


	nbins = 10
	fig, axs = plt.subplots(1, 3)
	axs[0].hist(results[gene]['flight']['normalized'], bins=nbins)
	axs[0].set_title('flight')
	axs[1].hist(results[gene]['ground']['normalized'], bins=nbins)
	axs[1].set_title('ground')
	axs[2].hist(results[gene]['nonflight']['normalized'], bins=nbins)
	axs[2].set_title('non-flight')
	fig.suptitle(gene)
	plt.savefig(gene + '.png', dpi=300)

	colors = ['red', 'blue', 'orange']
	combo_data = [results[gene]['flight']['normalized'], results[gene]['ground']['normalized'], results[gene]['nonflight']['normalized']]
	plt.hist(x=combo_data, bins=nbins, color=colors, density=True, stacked=True)
	plt.legend(['flight', 'ground', 'non-flight'])
	plt.title(gene)
	plt.savefig(gene + '.png', dpi=300)


	# qqplot
	'''stats.probplot(expr_per_gene_per_group[gene]['flight'])
	stats.probplot(expr_per_gene_per_group[gene]['ground'])
	stats.probplot(nonflight)
	pylab.show()'''
	print("gene = ", gene,
		  "\n\tflt mean = ", results[gene]['flight']['mean'],
		  "\n\tflt variance = ", results[gene]['flight']['variance'],
		  "\n\tbasal mean = ", results[gene]['basal']['mean'],
		  "\n\tbasal variance = ", results[gene]['basal']['variance'],
		  "\n\tground mean = ", results[gene]['ground']['mean'],
		  "\n\tground variance = ", results[gene]['ground']['variance'],
		  "\n\tvivarium mean = ", results[gene]['vivarium']['mean'],
		  "\n\tvivarium variance = ", results[gene]['vivarium']['variance'],
		  "\n\tnonflight mean = ", results[gene]['nonflight']['mean'],
		  "\n\tnonflight variance = ", results[gene]['nonflight']['variance'],
		  "\n\tflt_v_basal = ", results[gene]['flight_vs_basal'],
		  "\n\tflt_v_vivarium = ", results[gene]['flight_vs_vivarium'],
		  "\n\tflt_v_ground = ", results[gene]['flight_vs_ground'],
		  "\n\tflt_v_nonflight= ", results[gene]['flight_vs_nonflight'])

#import json
#print(json.dumps(results, sort_keys=False, indent=4))

