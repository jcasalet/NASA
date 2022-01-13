from sklearn.cluster import KMeans
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

expr_df=pd.read_csv('/Users/jcasaletto/Desktop/RESEARCH/NASA/BENCHMARK/DATA/expr_reduced_100_100.csv',
                    sep=',', header=0)

info_df=pd.read_csv('/Users/jcasaletto/Desktop/RESEARCH/NASA/BENCHMARK/DATA/all_metadata_Proj2.csv',
                    sep=',', header=0)

print("df shape = ", expr_df.shape)
# find max value
index=expr_df[expr_df['gene'] == expr_df.max()['gene']].index[0]
expr_df = expr_df.drop([expr_df.index[index]])

X=expr_df.drop(columns=['gene']).to_numpy()
print("X shape after drop = ", X.shape)

# find outliers
from scipy import stats
Z=np.abs(stats.zscore(X))
print(Z.sum())

km = KMeans(
    n_clusters=2, init='random',
    n_init=10, max_iter=300,
    tol=1e-04, random_state=0
)
y_km = km.fit_predict(X.T)
# plot the 3 clusters
plt.scatter(
    X.T[y_km == 0, 0], X.T[y_km == 0, 1],
    s=50, c='lightgreen',
    marker='s', edgecolor='black',
    label='cluster 1'
)

plt.scatter(
    X.T[y_km == 1, 0], X.T[y_km == 1, 1],
    s=5, c='orange',
    marker='o', edgecolor='black',
    label='cluster 2'
)

plt.show()