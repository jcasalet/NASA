from sklearn.decomposition import FastICA
import numpy as np
from matplotlib import pyplot as plt
from math import pi
from scipy.stats import ks_2samp
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

def create_ellipse(center_x, center_y, major_axis, minor_axis):
    t = np.linspace(0, 2 * pi, 500)
    x = [[center_x + major_axis * np.cos(i)] for i in t]
    y = [[center_y + minor_axis * np.sin(i)] for i in t]
    return np.hstack((x, y))


def add_noise(mu, sigma, X):
    for row in X:
        x = row[0]
        x_noise = np.random.normal(mu, sigma, 1)[0]
        row[0] = x + x_noise
        y = row[1]
        y_noise = np.random.normal(mu, sigma, 1)[0]
        row[1] = y + y_noise
    return X


def stretch(alpha, X):
    S = [[alpha, 0], [0, alpha]]
    return np.dot(X, S)


def rotate(theta, X):
    Y = np.zeros(X.shape)
    R = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    for i in range(len(X)):
        row = X[i]
        Y[i] = np.dot(R, row)
    return Y


X1 = create_ellipse(center_x=3, center_y=-1, major_axis=10, minor_axis=0.125)
X1 = rotate(90, X1)
X1 = add_noise(mu=0, sigma=2.0, X=X1)
X1 = stretch(2, X1)

X2 = create_ellipse(center_x=-3, center_y=1, major_axis=10, minor_axis=0.125)
X2 = rotate(180, X2)
X2 = add_noise(mu=0, sigma=1.5, X=X2)
X2 = stretch(2, X2)

X3 = create_ellipse(center_x=0, center_y=0, major_axis=10, minor_axis=0.125)
X3 = rotate(0, X3)
X3 = add_noise(mu=0, sigma=1, X=X3)
X3 = stretch(2, X3)

X1_X2 = np.concatenate((X1, X2), axis=0)
X1_X2_X3 = np.concatenate((X1_X2, X3), axis=0)

# find 3 independent components
# ica = FastICA(n_components=3, max_iter=1000)
ica = FastICA(max_iter=1000, n_components=2, random_state=0)

# run ICA on X
S = ica.fit_transform(X1_X2_X3)
S_inverse = ica.inverse_transform(S, copy=True)

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5)
fig.suptitle('original, combined, and reconstructed distributions')

for row in X1:
    ax1.scatter(row[0], row[1], color='yellow', marker='o', s=1)

for row in X2:
    ax1.scatter(row[0], row[1], color='red', marker='o', s=1)

for row in X3:
    ax1.scatter(row[0], row[1], color='green', marker='o', s=1)

for row in X1_X2_X3:
    ax2.scatter(row[0], row[1], color='orange', marker='o', s=1)

for row in S:
    ax3.scatter(row[0], row[1], color='black', marker='o', s=1)

for row in S_inverse:
    ax4.scatter(row[0], row[1], color='orange', marker='o', s=1)

S_dot_A = np.dot(S, ica.mixing_)
for row in S_dot_A:
    ax5.scatter(row[0], row[1], color='blue', marker='o', s=1)

#plt.savefig('./ica_with_3_dists.png', dpi=300)

plt.show()

plt.close()

print('params = ' + str(ica.get_params()))

print('S = ' + str(S))
print('components = ' + str(ica.components_))
print('mixing = ' + str(ica.mixing_))

# linearize X1_X2_X3
xlist = list()
slist = list()
rows = X1_X2_X3.shape[0]
cols = X1_X2_X3.shape[1]
for i in range(rows):
    for j in range(cols):
        xlist.append(X1_X2_X3[i][j])
        slist.append(S_dot_A[i][j])

print('2 sample test = ' + str(ks_2samp(xlist, slist, mode='asymp')))

# PCA
pca = PCA(n_components=2)
pca.fit(X1_X2_X3)
print(pca.explained_variance_ratio_)

# pearson
p_r = pearsonr(xlist, slist)
print(p_r)

