import pandas as pd
import numpy as np
import sys

fakeDFFile = sys.argv[1]
realDFFile = sys.argv[2]

fake_df = pd.read_csv(fakeDFFile, sep=',', header=0)
real_df = pd.read_csv(realDFFile, sep=',', header=0)

real_df = real_df.drop(columns=['gene'])
real_std = real_df.std()
real_mean = real_df.mean()

genes = fake_df['gene']

fake_df = fake_df.drop(columns=['gene'])
theMin = min(fake_df.min())
print(theMin)
if theMin < 0:
    theMin = -1 * theMin
    fake_df = fake_df.add(theMin)

fake_df = np.power(fake_df, 10)
fake_df = fake_df * real_std + real_mean

fake_df.insert(0, 'gene', genes)
fake_df.to_csv(fakeDFFile.split('.csv')[0] + '_shifted.csv', sep=',', index=False)

