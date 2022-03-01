import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file name', default=None, required=True)
    parser.add_argument('-m', '--meta', help='metadata file name', default=None, required=True)
    parser.add_argument('-o', '--out', help='output file name', default=None, required=True)
    parser.add_argument('-ce', '--crisp_env', help='colon-sep list of vars for crisp env', default=None, required=True)
    parser.add_argument('-k', '--key', help='key to index sample', default=None, required=True)
    return parser.parse_args()

def create_dict(keyName, sampleList, metaDF):
    myDict = dict()
    for sample in sampleList:
        value = metaDF[metaDF['sample'] == sample][keyName].iloc[0]
        myDict[sample] = value
    return myDict

def main():

    args = parse_args()
    # read in expr
    df = pd.read_csv(args.expr, header=0, sep=',')
    # read in meta
    metaDF = pd.read_csv(args.meta, header=0, sep=',')

    genes=list(df['gene'])
    df=df.drop(columns=['gene'])
    df = df.transpose()
    df.columns = genes
    df.reset_index(inplace=True)
    df = df.rename(columns={"index": "sample"})

    # replace "." with "_" in sample names
    #df['sample'] = df['sample'].str.replace('.', '-')

    envVar = args.crisp_env
    keyVar = args.key
    envVarList = []
    envVarArray = envVar.split(':')
    for e in envVarArray:
        envVarList.append(e)

    print(envVarList)
    # create env dict
    env_dict = dict()
    for i in range(len(metaDF)):
        sample = metaDF.iloc[i][keyVar]
        envString = ""
        for v in envVarList:
            if len(envString) == 0:
                envString = metaDF.iloc[i][v]
            else:
                envString = envString + ':' + metaDF.iloc[i][v]
        env_dict[sample] = envString
    df['env'] = df[keyVar].map(env_dict)

    # create other dicts
    for key in metaDF.columns:
        print('processing ' + key)
        if key == 'sample':
            continue
        theDict = create_dict(key, list(df['sample']), metaDF)
        # join dict to dfs
        df[key] = df['sample'].map(theDict)



    # save final uber table to pickle file
    df.to_pickle(args.out + '.pkl')
    df.to_csv(args.out + '.csv', sep=',', index=None)

if __name__ == "__main__":
    main()