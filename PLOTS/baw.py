import pandas as pd
import matplotlib.pyplot as plt
import argparse
import scipy.stats as stats



parser = argparse.ArgumentParser()
parser.add_argument('-m', '--meta', help='metadata file', default=None, required=True)
parser.add_argument('-f', '--field', help='metadata field', default=None, required=True)
args = parser.parse_args()

df=pd.read_csv(args.meta, header=0, sep=',')

field = args.field

fieldValues = set(df[field])


value_dict=dict()
results = dict()


flight = str(field) + '_flight'
nonflight= str(field) + '_nonflight'
results[field] = dict()
value_dict[flight] = list(df[df['Factor Value: Spaceflight'] == 'Space Flight'][field]) 
value_dict[nonflight] = list(df[df['Factor Value: Spaceflight'] == 'Ground Control'][field]) 

if len(value_dict[flight]) != 0 and len(value_dict[nonflight]) != 0:
    results[field]['t-test'] = float('%.5f' % (stats.ttest_ind(value_dict[flight], value_dict[nonflight]).pvalue))
    results[field]['wilcoxon'] = float('%.5f' % (stats.ranksums(value_dict[flight], value_dict[nonflight]).pvalue))
    results[field]['ks-test'] = float('%.5f' % (stats.kstest(value_dict[flight], value_dict[nonflight]).pvalue))


print(results)
fig,ax = plt.subplots()
ax.boxplot(value_dict.values())
ax.set_xticklabels(value_dict.keys())
#fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_' + field + '.png')

