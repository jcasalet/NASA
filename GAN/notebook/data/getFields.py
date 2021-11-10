import pandas as pd
import numpy as np

df = pd.read_csv('liver_meta_data.csv') 

# sample,Study,Group,oro,Dissection,strain,duration,library,GLDS

vocab_dicts = { 
'Study': {None: 0, 'CASIS': 1, 'RR-3': 2, 'RR1': 3},
'Group': {None: 0, 'Basal': 1, 'Flight': 2, 'Ground': 3, 'Vivarium': 4},
'Dissection': {None: 0, 'Frozen carcass': 1, 'Immediate': 2},
'strain': {None: 0, 'BALB/cT': 1, 'C57BL/6T': 2},
'library': {None: 0, 'polyA':1, 'ribo-depleted': 2}}

columns=['Study', 'Group', 'Dissection', 'strain', 'library']

encoded_conditions = []
for i in range(len(df)):
	myrow = []
	for col in columns:
		myval = vocab_dicts[col][df.iloc[i][col]]
		myrow.append(myval)
	encoded_conditions.append(np.array(myrow))



###################

expr_df = pd.read_csv('liver_with_flight.tsv', sep='\t', header=0)
sample_names = expr_df['sample'] 
gene_names = expr_df.columns 

colsToDrop = ['sample', 'env', 'oro_thresh', 'flight']
expr_df = expr_df.drop(columns=colsToDrop)

expression_values = []
for i in range(len(expr_df)):
	expr = list()
	for col in expr_df.columns:
		expr.append(expr_df[col])
	expression_values.append(expr)
expression_values = np.array(expression_values, dtype=np.float64).T

#print(expression_values)
print(gene_names)
