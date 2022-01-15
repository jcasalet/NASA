import pandas as pd
import numpy as np
import sys

dfFile = sys.argv[1]

df = pd.read_csv(dfFile, sep=',', header=0)

genes = df['gene']

df = df.drop(columns=['gene'])
theMin = min(df.min())
print(theMin)
if theMin < 0:
    theMin = -1 * theMin
    df = df.add(theMin)

df.insert(0, 'gene', genes)
df.to_csv(dfFile.split('.csv')[0] + '_min.csv', sep=',', index=False)

