import numpy as np
import sys
import pandas as pd


def pearson_correlation(x, y, kappa):
    """
    Computes similarity measure between each pair of genes in the bipartite graph x <-> y
    :param x: Gene matrix 1. Shape=(nb_samples, nb_genes_1)
    :param y: Gene matrix 2. Shape=(nb_samples, nb_genes_2)
    :return: Matrix with shape (nb_genes_1, nb_genes_2) containing the similarity coefficients
    """

    def standardize(a, kappa=1):
        a_off = np.mean(a, axis=0)
        a_std = np.std(a, axis=0)
        if not np.any(a_std):
            return np.zeros(a.shape)
        return (a - a_off) / (a_std * kappa)

    assert x.shape[0] == y.shape[0]
    x_ = standardize(x, kappa)
    y_ = standardize(y, kappa)
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


def correlations_list(x, y, corr_fn=pearson_correlation, kappa=1):
    """
    Generates correlation list between all pairs of genes in the bipartite graph x <-> y
    :param x: Gene matrix 1. Shape=(nb_samples, nb_genes_1)
    :param y: Gene matrix 2. Shape=(nb_samples, nb_genes_2)
    :param corr_fn: correlation function taking x and y as inputs
    """
    #corr = corr_fn(x, y, kappa)
    #return upper_diag_list(corr)
    return upper_diag_list(pearson_correlation(x, y, kappa))


def gamma_coef(x, y, kappa=1):
    """
    Compute gamma coefficients for two given expression matrices
    :param x: matrix of gene expressions. Shape=(nb_samples_1, nb_genes)
    :param y: matrix of gene expressions. Shape=(nb_samples_2, nb_genes)
    :return: Gamma(D^X, D^Z)
    """
    dists_x = 1 - correlations_list(x, x, kappa)
    dists_y = 1 - correlations_list(y, y, kappa)
    gamma_dx_dy = pearson_correlation(dists_x, dists_y, kappa)
    return gamma_dx_dy

def main():
    df1File = sys.argv[1]
    df2File = sys.argv[2]
    num_samples = int(sys.argv[3])

    df1 = pd.read_csv(df1File, sep=',', header=0).drop(columns=['gene'])
    df2 = pd.read_csv(df2File, sep=',', header=0).drop(columns=['gene'])

    print('df1 shape: ', str(df1.shape))
    print('df2 shape: ', str(df2.shape))

    df1_np = df1.to_numpy().T[0:num_samples].T
    # Log-transform data
    x = np.log(1 + df1_np)
    x = np.float32(x)
    # standardize expression data
    x = (x - x.mean()) / x.std()

    df2_np = df2.to_numpy()

    print('gamma = ', str(gamma_coef(x, df2_np)))

if __name__ == "__main__":
    main()