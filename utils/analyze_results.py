import pandas as pd
from statistics import median, mean
import sys

results_json_file = sys.argv[1]

results=pd.read_json(results_json_file)['results']

for i in range(len(results)):
	if 'method' in results.iloc[i].keys():
		print(results.iloc[i]['method'])
		print('min coeff = ', min(results.iloc[i]['coefficients']))
		print('max coeff = ', max(results.iloc[i]['coefficients']))
		print('mean coeff = ', mean(results.iloc[i]['coefficients']))

