from sklearn.impute import KNNImputer
from sklearn.impute import SimpleImputer
import sys
import pandas as pd
import numpy as np

inputFile = sys.argv[1]
fieldSep = sys.argv[2]
outputFile = sys.argv[3]

inputDF = pd.read_csv(inputFile, header=0, sep=fieldSep)

# remove any columns with all NaNs
allNanCols = inputDF.loc[:, inputDF.isnull().all()].columns.tolist()
print(allNanCols)
nonNullDF = inputDF.drop(columns=allNanCols)
nonNullColumns = nonNullDF.columns

# convert to numpy array
inputArray = nonNullDF.to_numpy()

# impute missing values
imputer = KNNImputer(n_neighbors=3, weights='uniform')
#imputer = SimpleImputer(missing_values=np.NaN, strategy='median', verbose=100)
#imputer = SimpleImputer(missing_values=np.NaN, strategy='constant', verbose=100)
imputedArray = imputer.fit_transform(inputArray)

# add cols back to get a df
imputedDF = pd.DataFrame(data=imputedArray,  columns=nonNullColumns)

# write df to file
imputedDF.to_csv(outputFile, index=False, sep=',')




