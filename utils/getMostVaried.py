import operator
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ie', '--input_expr', help='input expr file', default=None)
    parser.add_argument('-n', '--num', help='number of genes', default=0)
    parser.add_argument('-oe', '--output_expr', help='output expr file', default=None)
    return parser.parse_args()

# df is genes X samples
# calculate var, sort cols into n highest vars, drop shape[1]-n cols

def main():
    options = parse_args()
    n=int(options.num)
    exprFile=options.input_expr
    outputFile=options.output_expr
    df=pd.read_csv(exprFile, sep=',', header=0, index_col='gene')

    if n == 0:
        return df, None

    sdList = df.std(axis=1)
    sdDict = {k: v for v, k in enumerate(sdList)}
    sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=True)
    topN = sdDictSorted[0:n]
    indices = [x[1] for x in topN]
    #slicedDF = df[:,indices]
    slicedDF = df.iloc[indices]
    slicedDF.to_csv(outputFile, sep=',')

if __name__ == "__main__":
    main()