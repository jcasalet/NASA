iimport pandas as pd
from statistics import median, mean
import sys
import matplotlib.pyplot as plt

results_json_file = sys.argv[1]

results=pd.read_json(results_json_file)['results']

for i in range(len(results)):
        if 'method' in results.iloc[i].keys():
                method=results.iloc[i]['method']
                coefs=list(results.iloc[i]['coefficients'])
                print(method)
                print('min coeff = ', min(coefs))
                print('max coeff = ', max(coefs))
                print('mean coeff = ', mean(coefs))
                plt.hist(coefs)
                plt.savefig(method + '_hist.png')
