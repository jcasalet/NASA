import pandas as pd
import json
import sys

df=pd.read_csv(sys.argv[1])
genes=dict()

for i in range(len(df.columns)-1):
	genes[i]=list(df.sort_values(by=str(i), ascending=False)['gene'])[0:20]

with open('genes.json', 'w') as json_file:
	json.dump(genes, json_file)
json_file.close()
