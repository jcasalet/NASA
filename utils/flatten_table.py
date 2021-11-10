'''
SID,RADTYPE,DOSE,COHORT,MONTH,HR,LVAW_d,LVAW_s,LVID_d,LVID_s,LVPW_d,LVPW_s,LVmass,D_s,D_d,V_s,V_d,SV,EF,FS,CO,MVA,MVE,MVEA
11,2,0.1,2,0,507.154,0.934,1.408,3.24,1.668,0.89,1.537,100.993,1.808,3.236,9.88,42.126,32.246,76.589,44.137,16.353,,,
11,2,0.1,2,3,462.242,0.979,1.292,3.794,2.742,0.91,1.284,135.789,2.604,3.799,24.706,61.915,37.209,60.099,31.451,17.2,,,
11,2,0.1,2,7,390.476,0.784,1.101,4.302,3.303,0.756,1.032,126.086,3.248,4.313,42.465,83.684,41.218,49.254,24.709,16.099,,,
11,2,0.1,2,9,385.131,0.8,1.198,4.282,3.218,0.76,1.048,127.268,3.135,4.329,39.022,84.457,45.434,53.861,27.614,17.499,,,
12,2,0.1,2,0,457.425,1.121,1.679,4.052,2.627,0.788,1.184,153.009,2.596,3.992,24.525,69.682,45.156,64.802,34.974,20.655,,,
12,2,0.1,2,3,352.971,0.865,1.15,4.148,3.4,0.731,0.934,124.456,3.394,4.187,47.301,78.054,30.753,39.469,18.954,10.85,,,
12,2,0.1,2,7,412.934,0.999,1.235,4.408,3.583,0.756,0.979,156.293,3.487,4.451,50.423,90.118,39.695,43.994,21.643,16.399,,,
12,2,0.1,2,9,404.367,0.825,1.06,4.31,3.343,0.776,1.04,133.144,3.313,4.354,44.841,85.997,41.156,47.926,23.915,16.691,,,
13,2,0.1,2,0,488.978,1.028,1.365,3.599,2.462,0.808,1.296,120.01,2.357,3.576,19.356,53.669,34.314,63.958,34.08,16.754,,,
'''
import pandas as pd
import sys
import numpy as np

inputTable = sys.argv[1]
fieldDelimiter = sys.argv[2]
identifierField = sys.argv[3]
pivotColumn = sys.argv[4]

# converters={'SID': str,'MONTH':str}
inputDF = pd.read_csv(inputTable, sep=fieldDelimiter, header=0, skip_blank_lines=True).sort_values('SID')
inputDF.dropna(how="all", inplace=True)
inputDF = inputDF.astype({'SID': 'int', 'MONTH': 'str'})

columnPrefixes = inputDF.columns.drop([identifierField, pivotColumn])

flattenedDict = dict()

total = len(inputDF)
i = 0
while i < total:
    id = inputDF.iloc[i][identifierField]
    if id == '' or id is np.nan:
        continue
    idSubset = inputDF[inputDF[identifierField] == id]
    flattenedDict[id] = dict()
    numRows = len(idSubset)
    for j in range(numRows):
        pc = idSubset.iloc[j][pivotColumn]
        if not '0.' in pc:
            pc = str(int(float(pc)))
        for col in columnPrefixes:
            flattenedDict[id][str(col) + '_' + pc] = idSubset.iloc[j][col]
    i += numRows

flatDF = pd.DataFrame.from_dict(flattenedDict).T.reset_index()
flatDF=flatDF.fillna('NaN')
print('shape before: ' + str(flatDF.shape))
flatDF = flatDF.dropna(axis=1, how='all')
print('shape after: ' + str(flatDF.shape))
flatDF.to_csv(inputTable.replace('.csv', '.tsv'), sep='\t', index=False)
print('cols = ' + str(flatDF.columns))

