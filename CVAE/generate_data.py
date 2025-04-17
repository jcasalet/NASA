import torch
import requests
import numpy as np
import torch.nn.functional as F
from torch import nn, optim
import src.utils as utils
import scanpy as sc
from scipy import sparse
import matplotlib.pyplot as plt
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
import plotly.express as px
import plotly.io as pio
import pandas as pd
import src.vanilla_cvae as vanilla_cvae
import src.gma_cvae as gma_cvae
import seaborn as sns
from anndata import AnnData, concat
from scipy.stats import pearsonr

# Read in integrated data
adata = sc.read_h5ad('data/corrected_data.h5ad')
# adata.obs = adata.obs.drop(columns=['dataset'])
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]
print(adata)

np.random.seed(20)
random_indices = np.random.choice(adata.var.shape[0], size=10, replace=False)
genes_to_plot = adata.var_names[random_indices]
indices = [adata.var_names.get_loc(gene) for gene in genes_to_plot]

# GPU enable or disable
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Decide which model to use: 'vanilla' for vanilla CVAE
model_for_gen = 'gma'

# Vanilla CVAE
if model_for_gen == 'vanilla':
    model = vanilla_cvae.Vanilla_CVAE(n_genes=2000, n_labels=5, latent_size=64, beta=0.01, lr=0.001, wd=0.1, device=device)
    # van_cvae.train_net(train_loader = train_dataloader, test_loader = val_dataloader, n_epochs = 30)
    model.load_state_dict(torch.load('trained_models/trained_vanilla_cvae.pt', map_location=device))
# GMA CVAE
else:
    gmt_dict = utils.read_gmt('data/BrainGMTv2_MouseOrthologs.gmt', min_g=0, max_g=2000)
    gm_mask = utils.create_pathway_mask(adata.var.index.tolist(), gmt_dict, n_labels=5, add_missing=5, fully_connected=True)
    model = gma_cvae.GMA_CVAE(n_labels=5, pathway_mask=gm_mask, lr=0.01)
    # gma_cvae.train_net(train_loader = train_dataloader, test_loader = val_dataloader, n_epochs = 30)
    model.load_state_dict(torch.load('trained_models/trained_gma_cvae.pt', map_location=device))

condition_cols = ['Strain', 'Sex', 'Age at Launch', 'Duration', 'Flight']

# generate data
existing_data = adata.X.A.copy()
existing_data_612 = adata[adata.obs['dataset'] == '612']
existing_data_612 = existing_data_612.X.A
print(existing_data_612)

# modify conditions to desired values ['Strain' = 0 - 2, 'Sex' = 0 (female), 'Age at Launch' = int of weeks, 
# 'Duration' = int of days, 'Flight' = 0 for ground, 1 for flight]

# Label mapping: {'B6129SF2/J': 0 (612), 'BALB/c': 1 (352), 'C57BL/6NTac': 2 (613)}
conditions = torch.tensor([0, 0, 14, 28, 1], dtype=torch.float32)
generated_data_612 = model.generate(conditions, num_samples = 5000).numpy()
# generated_data_612.obs['Strain'] = 'B6129SF2/J'

conditions = torch.tensor([2, 0, 29, 53, 1], dtype=torch.float32)
generated_data_613 = model.generate(conditions, num_samples = 5000).numpy()
# generated_data_613.obs['Strain'] = 'C57BL/6NTac'

conditions = torch.tensor([1, 0, 12, 41, 1], dtype=torch.float32)
generated_data_352 = model.generate(conditions, num_samples = 5000).numpy()
# generated_data_352.obs['Strain'] = 'BALB/c'

generated_data = np.concatenate((generated_data_612, generated_data_613, generated_data_352))

real_means = np.mean(existing_data_612, axis = 0)
real_std = np.std(existing_data_612, axis=0)
generated_means = np.mean(generated_data_612, axis=0)
generated_std = np.std(generated_data_612, axis=0)

# Plot average gene expression
fig, ax = plt.subplots(figsize=(10, 8))

real_means_to_plot = real_means[indices]
real_std_to_plot = real_std[indices]
real_vars = np.var(existing_data[indices], axis=0)
generated_means_to_plot = generated_means[indices]
generated_std_to_plot = generated_std[indices]
generated_vars = np.var(generated_data[indices], axis=0)

# Scatter plot for real vs. generated data
ax.scatter(real_means_to_plot, generated_means_to_plot, alpha=0.7)

# Label each point with the gene name
for i, gene in enumerate(genes_to_plot):
    ax.text(real_means_to_plot[i], generated_means_to_plot[i], gene, fontsize=12, ha='right')

ax.set_title('Average Gene Expression: Real vs Generated')
ax.set_xlabel('Real')
ax.set_ylabel('Generated')
ax.plot([real_means_to_plot.min(), real_means_to_plot.max()], [real_means_to_plot.min(), real_means_to_plot.max()], 'k--')

plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(14, 8))

bar_width = 0.35
index = np.arange(len(genes_to_plot))

# Bar plots for real and generated data
bar1 = ax.bar(index, real_means_to_plot, bar_width, yerr=real_std_to_plot, capsize=5, label='Real')
bar2 = ax.bar(index + bar_width, generated_means_to_plot, bar_width, yerr=generated_std_to_plot, capsize=5, label='Generated')

# Adding labels
ax.set_xlabel('Genes')
ax.set_ylabel('Average Expression')
ax.set_title('Average Gene Expression: Real vs Generated')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(genes_to_plot, rotation=45)
ax.legend()

plt.tight_layout()
plt.show()

# Compute Pearson's R and R^2
mean_r, _ = pearsonr(real_means, generated_means)
var_r, _ = pearsonr(real_vars, generated_vars)
mean_r2 = mean_r ** 2
var_r2 = var_r ** 2

mean_std = np.std([real_means, generated_means], axis=0).mean()
var_std = np.std([real_vars, generated_vars], axis=0).mean()

fig, ax = plt.subplots(figsize=(10, 6))

# Data to plot
categories = ['Mean Expression', 'Variance']
r2_values = [mean_r2, var_r2]
error_bars = [mean_std, var_std]

# Create bar plot with error bars
bar_width = 0.3
bars = ax.bar(categories, r2_values, yerr=error_bars, capsize=5, width=bar_width, color=['blue', 'orange'])

# Adding labels
ax.set_xlabel('Metrics')
ax.set_ylabel('R^2')
ax.set_title('R^2 Values for Mean and Variance of Gene Expression', pad=20)
ax.set_ylim(0, 1)

# Adding the R^2 values on top of the bars
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, yval + 0.02, round(yval, 4), ha='center', va='bottom', fontsize=12)

# Adjust the layout to make plots closer together
plt.tight_layout(pad=2)
plt.show()
# filtered_adata = adata[adata.obs[column_name] == value_to_filter]

# unique_elements, counts = np.unique(generated_data, return_counts=True)

# Print the unique elements and their counts
# for element, count in zip(unique_elements, counts):
#     print(f"Element {element} occurs {count} times")

combined_data = np.vstack([existing_data, generated_data])
labels = np.array(['Existing'] * existing_data.shape[0] + ['Generated'] * generated_data.shape[0])

pca = PCA(n_components=2)
pca_result = pca.fit_transform(combined_data)

# Create a DataFrame for Plotly
pca_df = pd.DataFrame(pca_result, columns=['PCA Component 1', 'PCA Component 2'])

# Plotting with Plotly
fig = px.scatter(
    pca_df, x='PCA Component 1', y='PCA Component 2',
    color=labels,
    labels={'PCA Component 1': 'PCA Component 1', 'PCA Component 2': 'PCA Component 2'},
    title='PCA: Existing v. Generated'
)
fig.update_traces(marker=dict(size=2.5))
fig.show()

pca = PCA(n_components=3)
pca_result = pca.fit_transform(combined_data)

# Create a DataFrame for Plotly
pca_df = pd.DataFrame(pca_result, columns=['PCA Component 1', 'PCA Component 2', 'PCA Component 3'])

# Plotting with Plotly
fig = px.scatter_3d(
    pca_df, x='PCA Component 1', y='PCA Component 2', z = 'PCA Component 3',
    color=labels,
    labels={'PCA Component 1': 'PCA Component 1', 'PCA Component 2': 'PCA Component 2'},
    title='PCA: Existing v. Generated'
)
fig.update_traces(marker=dict(size=2))
fig.show()

umap_model = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1)
umap_result = umap_model.fit_transform(combined_data)

# Perform KMeans clustering
num_clusters = 6  # Define the number of clusters you want
kmeans = KMeans(n_clusters=num_clusters, random_state=0)
cluster_labels = kmeans.fit_predict(umap_result)

umap_df = pd.DataFrame(umap_result, columns=['UMAP Component 1', 'UMAP Component 2'])
umap_df['Cluster'] = cluster_labels
umap_df['Data Type'] = ['Existing'] * existing_data.shape[0] + ['Generated'] * generated_data.shape[0]

# Plotting with Plotly
fig = px.scatter(
    umap_df, x='UMAP Component 1', y='UMAP Component 2',
    color=labels,
    labels={'Data Type': 'Data Type'},
    title='UMAP: Color by Type'
)
fig.update_traces(marker=dict(size=3.5))
# Show the plot
fig.show()

