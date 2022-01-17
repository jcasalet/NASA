import tensorflow as tf
import numpy as np
import umap
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import sys
import pandas as pd


def predict(cc, nc, gen, z=None, training=False):
    """
    Make predictions
    :param cc: Categorical covariates
    :param nc: Numerical covariates
    :param gen: Generator model
    :param z: Latent input
    :param training: Whether training
    :return: Sampled data
    """
    nb_samples = cc.shape[0]
    if z is None:
        z_dim = gen.input[0].shape[-1]
        z = tf.random.normal([nb_samples, z_dim])
    out = gen([z, cc, nc], training=training)
    if not training:
        return out.numpy()
    return out

def split_train_test(x, train_rate=0.75, seed=0):
    """
    Split data into a train and a test sets
    :param train_rate: percentage of training samples
    :return: x_train, x_test
    """
    nb_samples = x.shape[0]
    split_point = int(train_rate * nb_samples)
    x_train = x[:split_point]
    x_test = x[split_point:]
    return x_train, x_test

def tsne_2d(data, **kwargs):
    """
    Transform data to 2d tSNE representation
    :param data: expression data. Shape=(dim1, dim2)
    :param kwargs: tSNE kwargs
    :return:
    """
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

def my_load():
    exprFile = sys.argv[1]
    metaFile = sys.argv[2]

    expr_df = pd.read_csv(exprFile, sep=',', header=0)
    meta_df = pd.read_csv(metaFile, sep=',', header=0)

    return expr_df, meta_df


# Load dataset
expr_df, info_df = my_load()
x = expr_df.values.T
symbols = expr_df.index.values
sampl_ids = expr_df.columns.values
tissues = info_df['libPrep'].values
datasets = info_df['dataset'].values

# Log-transform data
x = np.float32(x[1:])
x = np.log(1 + x)

# Process categorical metadata
cat_dicts = []
tissues_dict_inv = np.array(list(sorted(set(tissues))))
tissues_dict = {t: i for i, t in enumerate(tissues_dict_inv)}
tissues = np.vectorize(lambda t: tissues_dict[t])(tissues)
cat_dicts.append(tissues_dict_inv)
dataset_dict_inv = np.array(list(sorted(set(datasets))))
dataset_dict = {d: i for i, d in enumerate(dataset_dict_inv)}
datasets = np.vectorize(lambda t: dataset_dict[t])(datasets)
cat_dicts.append(dataset_dict_inv)
cat_covs = np.concatenate((tissues[:, None], datasets[:, None]), axis=-1)
cat_covs = np.int32(cat_covs)
print('Cat covs: ', cat_covs.shape)

# Process numerical metadata
num_covs = np.zeros((x.shape[0], 1), dtype=np.float32)
print('Num covs: ', num_covs.shape)

# Train/test split
np.random.seed(0)
idx = np.arange(x.shape[0])
np.random.shuffle(idx)
x = x[idx, :]
num_covs = num_covs[idx, :]
cat_covs = cat_covs[idx, :]

x_train, x_test = split_train_test(x)
num_covs_train, num_covs_test = split_train_test(num_covs)
cat_covs_train, cat_covs_test = split_train_test(cat_covs)

# Normalise data
x_mean = np.mean(x_train, axis=0)
x_std = np.std(x_train, axis=0)
# x_train = standardize(x_train, mean=x_mean, std=x_std)
# x_test = standardize(x_test, mean=x_mean, std=x_std)

gen = tf.keras.models.load_model('./gen_liver.h5')

x_gen = predict(cc=cat_covs_test, nc=num_covs_test, gen=gen)
x_gen = x_gen*x_std + x_mean
x_gen = np.clip(x_gen, 0, a_max=None)

x_combined = np.concatenate((x_test, x_gen))
categories = ['real'] * x_test.shape[0] + ['fake'] * x_gen.shape[0]
tissues_test = [tissues_dict_inv[tidx] for tidx in cat_covs_test[:, 0]]
tissues_combined = tissues_test + tissues_test

emb_2d = umap.UMAP().fit_transform(x_combined)

plt.figure(figsize=(10, 10))
plot_tsne_2d(emb_2d, labels=np.array(categories), s=4)
plt.title('UMAP real/synthetic')

plt.figure(figsize=(10, 10))
plot_tsne_2d(emb_2d, labels=np.array(tissues_combined), s=4)
plt.title('UMAP tissue type')

np.savetxt("train.csv", x_train, delimiter=",", header=",".join(symbols))
np.savetxt("test.csv", x_test, delimiter=",", header=",".join(symbols))
np.savetxt("gen.csv", x_gen, delimiter=",", header=",".join(symbols))
