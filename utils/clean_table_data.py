import pandas as pd
import sys

inputTable = sys.argv[1]
sep = sys.argv[2]
inputDF = pd.read_csv(inputTable, sep=sep, header=0, dtype=float)

inputArray = inputDF.to_numpy()

for i in range(len(inputArray)):
    for j in range(len(inputArray[0])):
        if not isinstance(inputArray[i][j], (int, float, complex)):
            print('i=' + str(i) + ' j=' + str(j) + ' value = ' + str(inputArray[i][j]))