import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--meta', help='input metadata file', default=None)
    parser.add_argument('-f', '--field', help='metadata field', default=None)
    return parser.parse_args()


def main():
    options = parse_args()
    metaFile = options.meta
    metaField = options.field
    df=pd.read_csv(metaFile, header=0, sep=',')

    fieldValues = list(set(df[metaField]))

    oro_dict = dict()
    for value in fieldValues:
        oro_dict[value] = dict()
        fig, ax = plt.subplots()
        for group in ['Basal', 'Vivarium', 'Ground', 'Flight']:
            oro_dict[value][group] = list(df[(df[metaField] == value) & (df['group'] == group)]['ORO Positivity (%)'])
            ax.boxplot(oro_dict[value][group])
            #ax.set_xticklabels(oro_dict[value][group])
            fig.set_size_inches(12, 4)
            #plt.savefig('ground_flight_bawplot_mmwithin.png')
            plt.show()

    print(oro_dict)

if __name__ == "__main__":
	main()
