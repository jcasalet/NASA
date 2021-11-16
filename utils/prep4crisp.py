import pandas as pd
import sys

inputFile=sys.argv[1]
sep='\t'

# removing leading comma in csv
# read in csv
#df=pd.read_csv('Normalized.CRISP.Liver.Symbol.Filtered.Subset.log2plus1.051721.csv', header=0, sep=',')
df=pd.read_csv(inputFile, header=0, sep=sep)


# set index and transpose
#df = df.transpose()

# re-index to set index as column name
df.reset_index(inplace=True)

# rename "index" column to "sample"
df=df.rename(columns={"index": "sample"})

# replace "." with "_" in sample names
df['sample'] = pd.to_numeric(df['sample'], downcast='integer')

# create dictionary with sample id as key and concatenated column strings as value
env_dict=dict()
for i in range(len(df)):
    key=df.iloc[i]['sample']
    value=str(df.iloc[i]['RADTYPE'])
    env_dict[key] = value

# join env dictionary to data frame
df['env'] = df['sample'].map(env_dict)

dose_dict=dict()
for i in range(len(df)):
    key=df.iloc[i]['sample']
    dose=int(df.iloc[i]['DOSE'])
    if dose == 0:
        value=0
    else:
        value=1
    dose_dict[key] = value

# join fliglht dictionary to data frame
df['dose_binary'] = df['sample'].map(dose_dict)

# save df to file
#df.to_pickle('cardio_env-RADTYPE_target-dose_bin.pkl')
df=df.set_index('sample')
df.to_csv('cardio_RADTYPE_dosebin.csv', sep=',', index=False)