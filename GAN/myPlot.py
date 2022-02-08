import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np


def pcaPlot(pca, df, info_df, variable, title, gen_dir):
    pcaDF = pd.DataFrame(data=pca.fit_transform(df), columns=['PC 1', 'PC 2'])
    pcaDF.index = info_df.index
    pcaDF = pd.concat([pcaDF, info_df[['condition']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['dataset']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['libPrep']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['mission']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['seqFacility']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['strain']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['gender']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['preservation']]], axis=1)


    sns.set(style="whitegrid", font_scale=1.1)
    fig, ax = plt.subplots(figsize=(5,5))

    ax = sns.scatterplot(x=pcaDF['PC 1'], y=pcaDF['PC 2'], hue=pcaDF[variable], s=100)

    ax.set_xlabel('PC 1 ' + '(' + str(round(pca.explained_variance_ratio_[0]*100, 1)) + '% variance)', fontsize=15)
    ax.set_ylabel('PC 2 ' + '(' + str(round(pca.explained_variance_ratio_[1]*100, 1)) + '% variance)', fontsize=15)
    ax.set_title(title, fontsize=20)
    #plt.show()
    if gen_dir is None:
        gen_dir = '.'
    plt.savefig(gen_dir + '/' + title, dpi=300)
    plt.close()

def tsne_2d(data, **kwargs):
    """
    Transform data to 2d tSNE representation
    :param data: expression data. Shape=(dim1, dim2)
    :param kwargs: tSNE kwargs
    :return:
    """
    from sklearn.manifold import TSNE
    print('... performing tSNE')
    tsne = TSNE(n_components=2, **kwargs)
    return tsne.fit_transform(data)

def plot_tsne_2d(data, labels, **kwargs):
    """
    Plots tSNE for the provided data, coloring the labels
    :param data: expression data. Shape=(dim1, dim2)
    :param labels: color labels. Shape=(dim1,)
    :param kwargs: tSNE kwargs
    :return: matplotlib axes
    """
    dim1, dim2 = data.shape

    # Prepare label dict and color map
    label_set = set(labels)
    label_dict = {k: v for k, v in enumerate(label_set)}

    # Perform tSNE
    if dim2 == 2:
        # print('plot_tsne_2d: Not performing tSNE. Shape of second dimension is 2')
        data_2d = data
    elif dim2 > 2:
        data_2d = tsne_2d(data, **kwargs)
    else:
        raise ValueError('Shape of second dimension is <2: {}'.format(dim2))

    # Plot scatterplot
    for k, v in label_dict.items():
        plt.scatter(data_2d[labels == v, 0], data_2d[labels == v, 1],
                    label=v)
    plt.legend()
    return plt.gca()

def myPlot(x, x_gen, info_df, output_dir):

    # tsne plots
    import umap.umap_ as umap
    print('x shape = ' + str(x.shape))
    print('x_gen shape = ' + str(x_gen.shape))
    x_combined = np.concatenate((x, x_gen))
    categories = ['real'] * x.shape[0] + ['fake'] * x_gen.shape[0]
    #tissues_test = [info_df[tidx] for tidx in info_df[:, 0]]
    #tissues_combined = tissues_test + tissues_test
    emb_2d = umap.UMAP().fit_transform(x_combined)
    plt.figure(figsize=(10, 10))
    plot_tsne_2d(emb_2d, labels=np.array(categories), s=4)
    plt.title('UMAP real/synthetic')
    #plt.show()
    plt.savefig(output_dir + '/umap_real_v_synthetic.png', dpi=300)

    pca = PCA(n_components=2)

    #x = standardize(x)

    pcaPlot(pca, x, info_df, 'condition', 'Condition_Real_Dataset_' + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'seqFacility', 'Sequencing_Facility_Real_Dataset_' + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'dataset', 'GLDS_Dataset_Real_Dataset_'  + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'libPrep', 'Library_Prep_Real_Dataset_'  + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'mission', 'Mission_Real_Dataset_'  + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'strain', 'Strain_Real_Dataset_'  + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'gender', 'Gender_Real_Dataset_' + 'n=' + str(x.shape[0]), output_dir)
    pcaPlot(pca, x, info_df, 'preservation', 'Preservation_Real_Dataset_' + 'n=' + str(x.shape[0]), output_dir)


    pcaPlot(pca, x_gen, info_df, 'condition', 'Condition_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'seqFacility', 'Sequencing_Facility_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'dataset', 'GLDS_Dataset_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'libPrep', 'Library_Prep_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'mission', 'Mission_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'strain', 'Strain_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'gender', 'Gender_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)
    pcaPlot(pca, x_gen, info_df, 'preservation', 'Preservation_Fake_Dataset_' + 'n=' + str(x_gen.shape[0]), output_dir)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-ie', '--input_expr', help='input expression data', default=None)
    parser.add_argument('-im', '--input_meta', help='input meta data', default=None)
    parser.add_argument('-ge', '--gen_expr', help='generated expr file', default=None)
    parser.add_argument('-od', '--output_dir', help='output dir', default=None, required=True)

    return parser.parse_args()


def main():

    options = parse_args()

    expr_df = pd.read_csv(options.input_expr, index_col=0)
    gen_df = pd.read_csv(options.gen_expr, index_col=0)
    info_df = pd.read_csv(options.input_meta, index_col=0)
    myPlot(expr_df, gen_df, info_df, options.output_dir)


if __name__ == "__main__":
    main()
