import pandas as pd
import sys

inputFileList = sys.argv[1]
outputFileName = sys.argv[2]
inputFiles = list()
inputFileArray = inputFileList.split(',')
for i in range(len(inputFileArray)):
    inputDF = pd.read_csv(inputFileArray[i].strip(), sep='\t', header=0, converters={'index': int})
    inputFiles.append(inputDF)

mondoDF = pd.concat(inputFiles, axis=0, ignore_index=False)
mondoDF=mondoDF.fillna('NaN')
mondoDF = mondoDF.dropna(axis=1, how='all')

mondoDF.to_csv(outputFileName, sep='\t', index=False)

