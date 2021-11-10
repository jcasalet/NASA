import tensorflow as tf
import datetime

from utils import *
from tf_utils import *
from collections import Counter
import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import operator
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from scipy.cluster.hierarchy import linkage, cophenet, dendrogram
import pickle
import random
from statsmodels.stats.multitest import multipletests
from matplotlib.collections import EllipseCollection
import scipy
import sys
import argparse

tfk = tf.keras
tfkl = tf.keras.layers
CHECKPOINTS_DIR = 'checkpoints/'


# ------------------
# RNASeqDB
# ------------------


# ---------------------
# DATA UTILITIES
# ---------------------
def standardize(x, mean=None, std=None):
    """
    Shape x: (nb_samples, nb_vars)
    """
    if mean is None:
        mean = np.mean(x, axis=0)
    if std is None:
        std = np.std(x, axis=0)
    return (x - mean) / std


def split_train_test(x, train_rate=0.75, seed=0):
    """
    Split data into a train and a test sets
    :param train_rate: percentage of training samples
    :return: x_train, x_test
    """
    random.seed(seed)
    nb_samples = x.shape[0]
    split_point = int(train_rate * nb_samples)
    x_train = x[:split_point]
    x_test = x[split_point:]
    return x_train, x_test


def split_train_test_v2(x, sampl_ids, train_rate=0.75):
    """
    Avoids patient leak between train/test set
    Split data into a train and a test sets
    :param train_rate: percentage of training samples
    :return: x_train, x_test
    """
    nb_samples = x.shape[0]
    sample_ids_rev = np.array([s[::-1] for s in sampl_ids])
    split_point = int(train_rate * nb_samples)
    idxs = np.argsort(sample_ids_rev)
    sample_ids_rev_sorted = sample_ids_rev[idxs]

    p = split_point
    while p == nb_samples and sample_ids_rev_sorted[p - 1] == sample_ids_rev_sorted[p - 2]:
        p += 1
    if p == nb_samples:
        raise Exception('Error: Cannot split samples into train and test sets')

    x_train = x[idxs[:p]]
    x_test = x[idxs[p:]]
    sample_ids_train = sampl_ids[idxs[:p]]
    sample_ids_test = sampl_ids[idxs[p:]]

    return x_train, x_test, sample_ids_train, sample_ids_test


def split_train_test_v3(sample_names, train_rate=0.75, seed=0):
    """
    Split data into a train and a test sets keeping replicates within the same set
    :param sample_names: list of sample names
    :param train_rate: percentage of training samples
    :param seed: random seed
    :return: lists of train and test sample indices
    """
    # Set random seed
    random.seed(seed)

    # Find replicate segments
    replicate_ranges = {}
    for i, name in enumerate(sample_names):
        repl_nb = int(name.split('_')[-1][1:])
        n = len(replicate_ranges)
        if repl_nb == 1:  # First replicate sample
            replicate_ranges[n] = {'start': i, 'end': i}
        else:
            replicate_ranges[n - 1]['end'] = i

    # Permute unique samples
    unique_sample_idxs = list(replicate_ranges.keys())
    random.shuffle(unique_sample_idxs)

    # Split data
    nb_unique = len(unique_sample_idxs)
    split_point = int(train_rate * nb_unique)
    unique_train = unique_sample_idxs[:split_point]

    # Recover replicates
    train_idxs = []
    test_idxs = []
    for i in range(nb_unique):
        for j in range(replicate_ranges[i]['start'], replicate_ranges[i]['end'] + 1):
            if i in unique_train:
                train_idxs.append(j)
            else:
                test_idxs.append(j)

    assert len(set(train_idxs + test_idxs)) == len(sample_names)
    return train_idxs, test_idxs


# ------------------
# E. coli M3D
# ------------------



# ---------------------
# CORRELATION UTILITIES
# ---------------------

def pearson_correlation(x, y):
    """
    Computes similarity measure between each pair of genes in the bipartite graph x <-> y
    :param x: Gene matrix 1. Shape=(nb_samples, nb_genes_1)
    :param y: Gene matrix 2. Shape=(nb_samples, nb_genes_2)
    :return: Matrix with shape (nb_genes_1, nb_genes_2) containing the similarity coefficients
    """

    def standardize(a):
        a_off = np.mean(a, axis=0)
        a_std = np.std(a, axis=0)
        return (a - a_off) / a_std

    assert x.shape[0] == y.shape[0]
    x_ = standardize(x)
    y_ = standardize(y)
    return np.dot(x_.T, y_) / x.shape[0]


def cosine_similarity(x, y):
    """
    Computes cosine similarity between vectors x and y
    :param x: Array of numbers. Shape=(n,)
    :param y: Array of numbers. Shape=(n,)
    :return: cosine similarity between vectors
    """
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


def upper_diag_list(m_):
    """
    Returns the condensed list of all the values in the upper-diagonal of m_
    :param m_: numpy array of float. Shape=(N, N)
    :return: list of values in the upper-diagonal of m_ (from top to bottom and from
             left to right). Shape=(N*(N-1)/2,)
    """
    m = np.triu(m_, k=1)  # upper-diagonal matrix
    tril = np.zeros_like(m_) + np.nan
    tril = np.tril(tril)
    m += tril
    m = np.ravel(m)
    return m[~np.isnan(m)]


def correlations_list(x, y, corr_fn=pearson_correlation):
    """
    Generates correlation list between all pairs of genes in the bipartite graph x <-> y
    :param x: Gene matrix 1. Shape=(nb_samples, nb_genes_1)
    :param y: Gene matrix 2. Shape=(nb_samples, nb_genes_2)
    :param corr_fn: correlation function taking x and y as inputs
    """
    corr = corr_fn(x, y)
    return upper_diag_list(corr)


def gamma_coef(x, y):
    """
    Compute gamma coefficients for two given expression matrices
    :param x: matrix of gene expressions. Shape=(nb_samples_1, nb_genes)
    :param y: matrix of gene expressions. Shape=(nb_samples_2, nb_genes)
    :return: Gamma(D^X, D^Z)
    """
    dists_x = 1 - correlations_list(x, x)
    dists_y = 1 - correlations_list(y, y)
    gamma_dx_dy = pearson_correlation(dists_x, dists_y)
    return gamma_dx_dy


# ---------------------
# CLUSTERING UTILITIES
# ---------------------

# ---------------------
# PLOTTING UTILITIES
# ---------------------

# %load adversarial-gene-expression/tf_utils.py



# ------------------
# LIMIT GPU USAGE
# ------------------

def limit_gpu(gpu_idx=0, mem=2 * 1024):
    """
    Limits gpu usage
    :param gpu_idx: Use this gpu
    :param mem: Maximum memory in bytes
    """
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        # Restrict TensorFlow to only allocate 1GB of memory on the first GPU
        try:
            # Use a single gpu
            tf.config.experimental.set_visible_devices(gpus[gpu_idx], 'GPU')

            # Limit memory
            tf.config.experimental.set_virtual_device_configuration(
                gpus[0],
                [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=mem)])  # 2 GB
            logical_gpus = tf.config.experimental.list_logical_devices('GPU')
            print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
        except RuntimeError as e:
            # Virtual devices must be set before GPUs have been initialized
            print(e)


# ------------------
# WGAN-GP
# ------------------

def make_generator(x_dim, vocab_sizes, nb_numeric, h_dims=None, z_dim=10):
    """
    Make generator
    :param x_dim: Number of genes
    :param vocab_sizes: List of ints. Size of vocabulary for each categorical covariate
    :param nb_numeric: Number of numerical covariates
    :param h_dims: Number of units for each hidden layer
    :param z_dim: Number of input units
    :return: generator
    """
    # Define inputs
    z = tfkl.Input((z_dim,))
    nb_categoric = len(vocab_sizes)
    cat = tfkl.Input((nb_categoric,), dtype=tf.int32)
    num = tfkl.Input((nb_numeric,), dtype=tf.float32)

    embed_cats = []
    total_emb_dim = 0

    for n, vs in enumerate(vocab_sizes):
        emb_dim = int(vs ** 0.5) + 1  # Rule of thumb
        c_emb = tfkl.Embedding(input_dim=vs,  # Vocabulary size
                               output_dim=emb_dim  # Embedding size
                               )(cat[:, n])
        embed_cats.append(c_emb)
        total_emb_dim += emb_dim
    if nb_categoric == 1:
        embeddings = embed_cats[0]
    else:
        embeddings = tfkl.Concatenate(axis=-1)(embed_cats)
    embeddings = tfkl.Concatenate(axis=-1)([num, embeddings])
    total_emb_dim += nb_numeric

    def make_generator_emb(x_dim, emb_dim, h_dims=None, z_dim=10):
        if h_dims is None:
            h_dims = [256, 256]

        z = tfkl.Input((z_dim,))
        t_emb = tfkl.Input((emb_dim,), dtype=tf.float32)
        h = tfkl.Concatenate(axis=-1)([z, t_emb])
        for d in h_dims:
            h = tfkl.Dense(d)(h)
            h = tfkl.ReLU()(h)
        h = tfkl.Dense(x_dim)(h)
        model = tfk.Model(inputs=[z, t_emb], outputs=h)
        return model

    gen_emb = make_generator_emb(x_dim=x_dim,
                                 emb_dim=total_emb_dim,
                                 h_dims=h_dims,
                                 z_dim=z_dim)
    model = tfk.Model(inputs=[z, cat, num], outputs=gen_emb([z, embeddings]))
    model.summary()
    return model


def make_discriminator(x_dim, vocab_sizes, nb_numeric, h_dims=None):
    """
    Make discriminator
    :param x_dim: Number of genes
    :param vocab_sizes: List of ints. Size of vocabulary for each categorical covariate
    :param nb_numeric: Number of numerical covariates
    :param h_dims: Number of units for each hidden layer
    :return: discriminator
    """
    if h_dims is None:
        h_dims = [256, 256]

    x = tfkl.Input((x_dim,))
    nb_categoric = len(vocab_sizes)
    cat = tfkl.Input((nb_categoric,), dtype=tf.int32)
    num = tfkl.Input((nb_numeric,), dtype=tf.float32)

    embed_cats = []

    for n, vs in enumerate(vocab_sizes):
        emb_dim = int(vs ** 0.5) + 1  # Rule of thumb
        c_emb = tfkl.Embedding(input_dim=vs,  # Vocabulary size
                               output_dim=emb_dim  # Embedding size
                               )(cat[:, n])
        embed_cats.append(c_emb)

    if nb_categoric == 1:
        embeddings = embed_cats[0]
    else:
        embeddings = tfkl.Concatenate(axis=-1)(embed_cats)
    h = tfkl.Concatenate(axis=-1)([x, num, embeddings])
    for d in h_dims:
        h = tfkl.Dense(d)(h)
        h = tfkl.ReLU()(h)
    h = tfkl.Dense(1)(h)
    model = tfk.Model(inputs=[x, cat, num], outputs=h)
    return model


def wasserstein_loss(y_true, y_pred):
    """
    Wasserstein loss
    """
    return tf.reduce_mean(y_true * y_pred)


def generator_loss(fake_output):
    """
    Generator loss
    """
    return wasserstein_loss(-tf.ones_like(fake_output), fake_output)


def gradient_penalty(f, real_output, fake_output):
    """
    Gradient penalty of WGAN-GP
    :param f: discriminator function without sample covariates as input
    :param real_output: real data
    :param fake_output: fake data
    :return: gradient penalty
    """
    alpha = tf.random.uniform([real_output.shape[0], 1], 0., 1.)
    diff = fake_output - real_output
    inter = real_output + (alpha * diff)
    with tf.GradientTape() as t:
        t.watch(inter)
        pred = f(inter)
    grad = t.gradient(pred, [inter])[0]
    slopes = tf.sqrt(tf.reduce_sum(tf.square(grad), axis=1))  # real_output
    gp = tf.reduce_mean((slopes - 1.) ** 2)
    return gp


def discriminator_loss(real_output, fake_output):
    """
    Critic loss
    """
    real_loss = wasserstein_loss(-tf.ones_like(real_output), real_output)
    fake_loss = wasserstein_loss(tf.ones_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss


# Notice the use of `tf.function`
# This annotation causes the function to be "compiled".
@tf.function
def train_disc(x, z, cc, nc, gen, disc, disc_opt, grad_penalty_weight=10, p_aug=0, norm_scale=0.5):
    """
    Train critic
    :param x: Batch of expression data
    :param z: Batch of latent variables
    :param cc: Batch of categorical covariates
    :param nc: Batch of numerical covariates
    :param gen: Generator
    :param disc: Critic
    :param disc_opt: Critic optimizer
    :param grad_penalty_weight: Weight for the gradient penalty
    :return: Critic loss
    """
    bs = z.shape[0]
    nb_genes = gen.output.shape[-1]
    augs = np.random.binomial(1, p_aug, bs)

    with tf.GradientTape() as disc_tape:
        # Generator forward pass
        x_gen = gen([z, cc, nc], training=False)

        # Perform augmentations
        x_gen = x_gen + augs[:, None] * np.random.normal(0, norm_scale, nb_genes)
        x = x + augs[:, None] * np.random.normal(0, norm_scale, (bs, nb_genes))

        # Forward pass on discriminator
        disc_out = disc([x_gen, cc, nc], training=True)
        disc_real = disc([x, cc, nc], training=True)

        # Compute losses
        disc_loss = discriminator_loss(disc_real, disc_out) \
                    + grad_penalty_weight * gradient_penalty(lambda x: disc([x, cc, nc], training=True), x, x_gen)

    disc_grad = disc_tape.gradient(disc_loss, disc.trainable_variables)
    disc_opt.apply_gradients(zip(disc_grad, disc.trainable_variables))

    return disc_loss


@tf.function
def train_gen(z, cc, nc, gen, disc, gen_opt, p_aug=0, norm_scale=1):
    """
    Train generator
    :param z: Batch of latent variables
    :param cc: Batch of categorical covariates
    :param nc: Batch of numerical covariates
    :param gen: Generator
    :param disc: Critic
    :param gen_opt: Generator optimiser
    :return: Generator loss
    """
    bs = z.shape[0]
    nb_genes = gen.output.shape[-1]
    augs = np.random.binomial(1, p_aug, bs)

    with tf.GradientTape() as gen_tape:
        # Generator forward pass
        x_gen = gen([z, cc, nc], training=True)

        # Perform augmentations
        x_gen = x_gen + augs[:, None] * np.random.normal(0, norm_scale, (bs, nb_genes))

        # Forward pass on discriminator
        disc_out = disc([x_gen, cc, nc], training=False)

        # Compute losses
        gen_loss = generator_loss(disc_out)

    gen_grad = gen_tape.gradient(gen_loss, gen.trainable_variables)
    gen_opt.apply_gradients(zip(gen_grad, gen.trainable_variables))

    return gen_loss


def train(dataset, cat_covs, num_covs, z_dim, epochs, batch_size, gen, disc, score_fn, save_fn,
          gen_opt=None, disc_opt=None, nb_critic=5, verbose=True, checkpoint_dir='./checkpoints/cpkt',
          log_dir='./logs/', patience=10, p_aug=0, norm_scale=0.5):
    """
    Train model
    :param dataset: Numpy matrix with data. Shape=(nb_samples, nb_genes)
    :param cat_covs: Categorical covariates. Shape=(nb_samples, nb_cat_covs)
    :param num_covs: Numerical covariates. Shape=(nb_samples, nb_num_covs)
    :param z_dim: Int. Latent dimension
    :param epochs: Number of training epochs
    :param batch_size: Batch size
    :param gen: Generator model
    :param disc: Critic model
    :param gen_opt: Generator optimiser
    :param disc_opt: Critic optimiser
    :param score_fn: Function that computes the score: Generator => score.
    :param save_fn:  Function that saves the model.
    :param nb_critic: Number of critic updates for each generator update
    :param verbose: Print details
    :param checkpoint_dir: Where to save checkpoints
    :param log_dir: Where to save logs
    :param patience: Number of epochs without improving after which the training is halted
    """
    # Optimizers
    if gen_opt is None:
        gen_opt = tfk.optimizers.RMSprop(5e-4)
    if disc_opt is None:
        disc_opt = tfk.optimizers.RMSprop(5e-4)

    # Set up logs and checkpoints
    current_time = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    gen_log_dir = log_dir + current_time + '/gen'
    disc_log_dir = log_dir + current_time + '/disc'
    gen_summary_writer = tf.summary.create_file_writer(gen_log_dir)
    disc_summary_writer = tf.summary.create_file_writer(disc_log_dir)

    checkpoint_prefix = os.path.join(checkpoint_dir, 'ckpt')
    checkpoint = tf.train.Checkpoint(generator_optimizer=gen_opt,
                                     discriminator_optimizer=disc_opt,
                                     generator=gen,
                                     discriminator=disc)

    gen_losses = tfk.metrics.Mean('gen_loss', dtype=tf.float32)
    disc_losses = tfk.metrics.Mean('disc_loss', dtype=tf.float32)
    best_score = -np.inf
    initial_patience = patience

    for epoch in range(epochs):
        for i in range(0, len(dataset), batch_size):
            x = dataset[i: i + batch_size, :]
            cc = cat_covs[i: i + batch_size, :]
            nc = num_covs[i: i + batch_size, :]

            # Train critic
            disc_loss = None
            for _ in range(nb_critic):
                z = tf.random.normal([x.shape[0], z_dim])
                disc_loss = train_disc(x, z, cc, nc, gen, disc, disc_opt, p_aug=p_aug, norm_scale=norm_scale)
            disc_losses(disc_loss)

            # Train generator
            z = tf.random.normal([x.shape[0], z_dim])
            gen_loss = train_gen(z, cc, nc, gen, disc, gen_opt, p_aug=p_aug, norm_scale=norm_scale)
            gen_losses(gen_loss)

        # Logs
        with disc_summary_writer.as_default():
            tf.summary.scalar('loss', disc_losses.result(), step=epoch)
        with gen_summary_writer.as_default():
            tf.summary.scalar('loss', gen_losses.result(), step=epoch)

        # Save the model
        if epoch % 5 == 0:
            checkpoint.save(file_prefix=checkpoint_prefix)

            score = score_fn(gen)
            if score > best_score:
                print('Saving model ...')
                save_fn()
                best_score = score
                patience = initial_patience
            else:
                patience -= 1

            if verbose:
                print('Score: {:.3f}'.format(score))

        if verbose:
            print('Epoch {}. Gen loss: {:.2f}. Disc loss: {:.2f}'.format(epoch + 1,
                                                                         gen_losses.result(),
                                                                         disc_losses.result()))
        gen_losses.reset_states()
        disc_losses.reset_states()

        if patience == 0:
            break


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



#### Data


def my_find_mostvaried(df, n):
    # df is genes X samples
    # calculate var, sort cols into n highest vars, drop shape[1]-n cols
    if n == 0:
        return df
    sdList = df.std(axis=0)
    sdDict = {k: v for v, k in enumerate(sdList)}
    sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=True) 
    topN = sdDictSorted[0:n]
    indices = [x[1] for x in topN]
    slicedDF = df[:,indices]
    return slicedDF

def my_prep_data(n, expr_df, info_df):
    # Load dataset
    #expr_df, info_df = rnaseqdb_load()
    x = expr_df.values.T
    x = my_find_mostvaried(x, n)
    #symbols = expr_df.index.levels[0].values
    symbols = expr_df.index.values
    sampl_ids = expr_df.columns.values
    #     tissues = info_df['TISSUE_GTEX'].values
    #     datasets = info_df['DATASET'].values
    # sample,study,dissection,strain,libprep
    tissues = info_df['condition'] # not renaming the tissues variable for convenience LMS
    datasets = info_df['dataset']
    lib = info_df['libPrep']
    mission = info_df['mission']
    seqfac = info_df['seqFacility']

    # Log-transform data
    x = np.log(1 + x)
    x = np.float32(x)

    # Process categorical metadata
    cat_dicts = [] # big dict to hold all categorical dicts
    def cat(var):
        '''Function to repeatedly process categorical metadata. Pass in a column ("var") from info_df as a variable.'''
        var_dict_inv = np.array(list(sorted(set(var))))
        var_dict = {t: i for i, t in enumerate(var_dict_inv)}
        var = np.vectorize(lambda t: var_dict[t])(var) # convert to integer
        cat_dicts.append(var_dict_inv) # add to big dict
        return var, var_dict_inv

    tissues, tissues_dict_inv = cat(tissues)
    datasets, datasets_dict_inv = cat(datasets)
    lib, lib_dict_inv = cat(lib)
    mission, mission_dict_inv = cat(mission)
    seqfac, seqfac_dict_inv = cat(seqfac)

    ## Final concatenation
    cat_covs = np.concatenate((tissues[:, None],
                               datasets[:, None],
                              lib[:, None],
                                mission[:, None],
                                seqfac[:, None]),
                              axis=-1)

    #print(cat_covs)
    cat_covs = np.int32(cat_covs) # make sure all are integers
    print('Cat covs: ', cat_covs.shape)
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
    x_train = standardize(x_train, mean=x_mean, std=x_std)
    x_test = standardize(x_test, mean=x_mean, std=x_std)

    return cat_dicts, cat_covs, cat_covs_test, cat_covs_train, num_covs, num_covs_test, num_covs_train, x, x_test, x_train


def my_train(CONFIG, cat_dicts, cat_covs, cat_covs_test, cat_covs_train, num_covs, num_covs_test, num_covs_train, x, x_test, x_train):
    # Train on GL liver...

    # Script that trains the model
    # %load adversarial-gene-expression/gtex_tcga_gan.py
    #import wandb

    #tf.config.run_functions_eagerly(True)  # LMS added to avoid variable creation error https://stackoverflow.com/questions/58352326/running-the-tensorflow-2-0-code-gives-valueerror-tf-function-decorated-functio

    MODELS_DIR = 'checkpoints/models/'

    # GPU limit
    limit_gpu(CONFIG['gpu'])

    # Define model
    vocab_sizes = [len(c) for c in cat_dicts]
    print('Vocab sizes: ', vocab_sizes)
    nb_numeric = num_covs.shape[-1]
    x_dim = x.shape[-1]
    gen = make_generator(x_dim, vocab_sizes, nb_numeric,
                         h_dims=[CONFIG['hdim']] * CONFIG['nb_layers'],
                         z_dim=CONFIG['latent_dim'])
    disc = make_discriminator(x_dim, vocab_sizes, nb_numeric,
                              h_dims=[CONFIG['hdim']] * CONFIG['nb_layers'])

    # Evaluation metrics
    def score_fn(x_test, cat_covs_test, num_covs_test):
        def _score(gen):
            x_gen = predict(cc=cat_covs_test, ## x_gen is an array of nans, throws downstream Assertion Error LMS
                            nc=num_covs_test,
                            gen=gen)

            gamma_dx_dz = gamma_coef(x_test, x_gen)
            return gamma_dx_dz
            #score = (x_test - x_gen) ** 2
            #return -np.mean(score)

        return _score

    # Function to save models
    def save_fn(models_dir=MODELS_DIR):
        gen.save(models_dir + 'gen_liver.h5')


    # Train model
    gen_opt = tfk.optimizers.RMSprop(CONFIG['lr'])
    disc_opt = tfk.optimizers.RMSprop(CONFIG['lr'])

    ## Come back to this LMS (do I need a project name?)
#     run = wandb.init(project='adversarial_gene_expr', config=CONFIG)
#     config = wandb.config
#     # wandb.run.name = '{}'.format(wandb.run.name)
#     wandb.run.save()

    train(dataset=x_train,
          cat_covs=cat_covs_train,
          num_covs=num_covs_train,
          z_dim=CONFIG['latent_dim'],
          batch_size=CONFIG['batch_size'],
          epochs=CONFIG['epochs'],
          nb_critic=CONFIG['nb_critic'],
          gen=gen,
          disc=disc,
          gen_opt=gen_opt,
          disc_opt=disc_opt,
          score_fn=score_fn(x_test, cat_covs_test, num_covs_test),
          save_fn=save_fn)

    # Evaluate data
    score = score_fn(x_test, cat_covs_test, num_covs_test)(gen)
    print('Gamma(Dx, Dz): {:.2f}'.format(score))


def pcaPlot(pca, df, info_df, variable, title):
    pcaDF = pd.DataFrame(data=pca.fit_transform(df), columns=['PC 1', 'PC 2'])
    pcaDF.index = info_df.index
    pcaDF = pd.concat([pcaDF, info_df[['condition']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['dataset']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['libPrep']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['mission']]], axis=1)
    pcaDF = pd.concat([pcaDF, info_df[['seqFacility']]], axis=1)
    sns.set(style="whitegrid", font_scale=1.1)
    fig, ax = plt.subplots(figsize=(5,5))

    ax = sns.scatterplot(x=pcaDF['PC 1'], y=pcaDF['PC 2'], hue=df[variable], s=100)

    ax.set_xlabel('PC 1 ' + '(' + str(round(pca.explained_variance_ratio_[0]*100, 1)) + '% variance)', fontsize=15)
    ax.set_ylabel('PC 2 ' + '(' + str(round(pca.explained_variance_ratio_[1]*100, 1)) + '% variance)', fontsize=15)
    ax.set_title(title, fontsize=20)
    #plt.show()
    plt.savefig('./' + title, dpi=300)
    plt.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--train', help='boolean train or not', default=True)
    parser.add_argument('-m', '--model', help='model file to use instead of training', default=None)
    parser.add_argument('-g', '--gpu', help='number of gpus', default=0)
    parser.add_argument('-e', '--epochs', help='epochs', default=100)
    parser.add_argument('-ld', '--latent_dim', help='number of latent dimensions', default=16)
    parser.add_argument('-bs', '--batch_size', help='batch size', default=8)
    parser.add_argument('-nl', '--nb_layers', help='number of layers', default=2)
    parser.add_argument('-hd', '--hdim', help='number of units per hidden layer ', default=256)
    parser.add_argument('-lr', '--lr', help='learning rate', default=5e-04)
    parser.add_argument('-nb', '--nb_critic', help='number of critic batches per gen batch', default=5)
    parser.add_argument('-ng', '--num_genes', help='number of genes with highest variance', default=0)
    return parser.parse_args() 
    
def main():
    options = parse_args()
    CONFIG = {'gpu': int(options.gpu), 'epochs': int(options.epochs), 'latent_dim': int(options.latent_dim),
              'batch_size': int(options.batch_size), 'nb_layers': int(options.nb_layers), 'hdim': int(options.hdim),
              'lr': float(options.lr), 'nb_critic': int(options.nb_critic)}
    data_dir = '.'
    genes_file = open("top-liver-genes.txt", "r")
    genes = genes_file. readlines()
    genes_list=list()
    for g in genes:
        genes_list.append(g.strip())
    expr_df = pd.read_csv(data_dir+'/Proj2_Normalized_Counts.csv', index_col=0)
    #df = pd.read_csv(data_dir+'/Normalized.CRISP.Liver.Symbol.Filtered.Subset.log2plus1.051721.csv', index_col=0).T
    #expr_df = df[df.columns.intersection(genes_list)].T
    info_df = pd.read_csv(data_dir+'/all_metadata_Proj2.csv', index_col=0)
    
    #standardize expression data
    expr_df = (expr_df-expr_df.mean())/expr_df.std()
    
    cat_dicts, cat_covs, cat_covs_test, cat_covs_train, num_covs, num_covs_test, num_covs_train, x, x_test, \
            x_train = my_prep_data(int(options.num_genes), expr_df, info_df)

    print('train option = ' + str(options.train))
    print('model file = ' + str(options.model))
    if eval(options.train):
        print('training!')
        my_train(CONFIG, cat_dicts, cat_covs, cat_covs_test, cat_covs_train, num_covs, num_covs_test, num_covs_train, x, \
             x_test, x_train)
        gen = tf.keras.models.load_model('checkpoints/models/gen_liver.h5') # this is the one I just trained
    else:
        print('not training!')
        gen = tf.keras.models.load_model(options.model)
    x_gen = predict(cc=cat_covs, nc=num_covs, gen=gen)
    print('x-gen shape = ' + str(x_gen.shape))
    print('x shape = ' + str(x.shape))
    pca = PCA(n_components=2)
    pcaPlot(pca, x, info_df, 'condition', 'Condition_Real_Dataset')
    pcaPlot(pca, x, info_df, 'seqFacility', 'Sequencing_Facility_Real_Dataset')
    pcaPlot(pca, x, info_df, 'dataset', 'GLDS_Dataset_Real_Dataset')
    pcaPlot(pca, x, info_df, 'libPrep', 'Library_Prep_Real_Dataset')
    pcaPlot(pca, x, info_df, 'mission', 'Mission_Real_Dataset')
    

    pcaPlot(pca, x_gen, info_df, 'condition', 'Condition_Fake_Dataset')
    pcaPlot(pca, x_gen, info_df, 'seqFacility', 'Sequencing_Facility_Fake_Dataset')
    pcaPlot(pca, x_gen, info_df, 'dataset', 'GLDS_Dataset_Fake_Dataset')
    pcaPlot(pca, x_gen, info_df, 'libPrep', 'Library_Prep_Fake_Dataset')
    pcaPlot(pca, x_gen, info_df, 'mission', 'Mission_Fake_Dataset')
                    

if __name__ == "__main__":
    main()
