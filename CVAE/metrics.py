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

van_model = vanilla_cvae.Vanilla_CVAE(n_genes=2000, n_labels=5, latent_size=64, beta=0.01, lr=0.001, wd=0.1, device=device)
# van_cvae.train_net(train_loader = train_dataloader, test_loader = val_dataloader, n_epochs = 30)
van_model.load_state_dict(torch.load('trained_models/trained_vanilla_cvae.pt', map_location=device))

# GMA CVAE
gmt_dict = utils.read_gmt('data/BrainGMTv2_MouseOrthologs.gmt', min_g=0, max_g=2000)
gm_mask = utils.create_pathway_mask(adata.var.index.tolist(), gmt_dict, n_labels=5, add_missing=5, fully_connected=True)
gma_model = gma_cvae.GMA_CVAE(n_labels=5, pathway_mask=gm_mask, lr=0.01)
# gma_cvae.train_net(train_loader = train_dataloader, test_loader = val_dataloader, n_epochs = 30)
gma_model.load_state_dict(torch.load('trained_models/trained_gma_cvae.pt', map_location=device))

condition_cols = ['Strain', 'Sex', 'Age at Launch', 'Duration', 'Flight']

# # generate data
existing_data = adata.X.A.copy()
existing_data_612 = adata[adata.obs['dataset'] == '612']
existing_data_612 = existing_data_612.X.A
print(existing_data_612)

# modify conditions to desired values ['Strain' = 0 - 2, 'Sex' = 0 (female), 'Age at Launch' = int of weeks, 
# 'Duration' = int of days, 'Flight' = 0 for ground, 1 for flight]

# Label mapping: {'B6129SF2/J': 0 (612), 'BALB/c': 1 (352), 'C57BL/6NTac': 2 (613)}
conditions = torch.tensor([0, 0, 14, 28, 1], dtype=torch.float32)
v_generated_data_612 = van_model.generate(conditions, num_samples = 5000).numpy()
g_generated_data_612 = gma_model.generate(conditions, num_samples = 5000).numpy()
# generated_data_612.obs['Strain'] = 'B6129SF2/J'

conditions = torch.tensor([2, 0, 29, 53, 1], dtype=torch.float32)
v_generated_data_613 = van_model.generate(conditions, num_samples = 5000).numpy()
g_generated_data_613 = gma_model.generate(conditions, num_samples = 5000).numpy()
# generated_data_613.obs['Strain'] = 'C57BL/6NTac'

conditions = torch.tensor([1, 0, 12, 41, 1], dtype=torch.float32)
v_generated_data_352 = van_model.generate(conditions, num_samples = 5000).numpy()
g_generated_data_352 = gma_model.generate(conditions, num_samples = 5000).numpy()
# generated_data_352.obs['Strain'] = 'BALB/c'

v_generated_data = np.concatenate((v_generated_data_612, v_generated_data_613, v_generated_data_352))
g_generated_data = np.concatenate((g_generated_data_612, g_generated_data_613, g_generated_data_352))


real_means = np.mean(existing_data, axis=0)
real_vars = np.var(existing_data, axis=0)

v_generated_means = np.mean(v_generated_data, axis=0)
v_generated_vars = np.var(v_generated_data, axis=0)

# Calculate means and variances for model 2 generated data
g_generated_means = np.mean(g_generated_data, axis=0)
g_generated_vars = np.var(g_generated_data, axis=0)

# Plot average gene expression
fig, ax = plt.subplots(figsize=(10, 8))

# Compute Pearson's R and R^2 for model 1
mean_r_v, _ = pearsonr(real_means, v_generated_means)
var_r_v, _ = pearsonr(real_vars, v_generated_vars)
mean_r2_v = mean_r_v ** 2
var_r2_v = var_r_v ** 2

# Compute Pearson's R and R^2 for model 2
mean_r_g, _ = pearsonr(real_means, g_generated_means)
var_r_g, _ = pearsonr(real_vars, g_generated_vars)
mean_r2_g = mean_r_g ** 2
var_r2_g = var_r_g ** 2

mean_std_v = np.std([real_means, v_generated_means], axis=0).mean()
var_std_v = np.std([real_vars, v_generated_vars], axis=0).mean()

# Calculate standard deviations for model 2
mean_std_g = np.std([real_means, g_generated_means], axis=0).mean()
var_std_g = np.std([real_vars, g_generated_vars], axis=0).mean()

# Plot R^2 values with error bars
fig, ax = plt.subplots(figsize=(12, 8))

# Data to plot
categories = ['Mean', 'Variance']
r2_values_v = [mean_r2_v, var_r2_v]
r2_values_g = [mean_r2_g, var_r2_g]
error_bars_v = [mean_std_v, var_std_v]
error_bars_g = [mean_std_g, var_std_g]

# Create bar plot with error bars
bar_width = 0.3
index = np.arange(len(categories))

bars_model1 = ax.bar(index, r2_values_v, yerr=error_bars_v, capsize=5, width=bar_width, color='blue', edgecolor='black', label='Vanilla CVAE')
bars_model2 = ax.bar(index + bar_width, r2_values_g, yerr=error_bars_g, capsize=5, width=bar_width, color='yellow', edgecolor='black', label='GMA CVAE')

# Adding labels
ax.set_xlabel('Metrics', fontsize=22)
ax.set_ylabel('R^2', fontsize=22)
ax.set_title('R^2 Values for Mean and Variance of Gene Expression\nBetween Real and Generated Cells', pad=20, fontsize=20)
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(categories, fontsize=22)
ax.set_ylim(0, 1.1)
ax.legend(fontsize=20)

# Adding the R^2 values on top of the bars
for bars in [bars_model1, bars_model2]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 0.02, round(yval, 4), ha='center', va='bottom', fontsize=18)

plt.show()

conditions = torch.tensor([0, 0, 14, 28, 1], dtype=torch.float32)
v_generated_data_flight = van_model.generate(conditions, num_samples = 5000).numpy()
g_generated_data_flight = gma_model.generate(conditions, num_samples = 5000).numpy()

conditions = torch.tensor([0, 0, 14, 28, 0], dtype=torch.float32)
v_generated_data_ground = van_model.generate(conditions, num_samples = 5000).numpy()
g_generated_data_ground = gma_model.generate(conditions, num_samples = 5000).numpy()

v_combined_data = np.vstack([existing_data, v_generated_data_612, v_generated_data_613, v_generated_data_352])
v_labels = np.array(['Existing'] * existing_data.shape[0] + ['Generated 612'] * v_generated_data_612.shape[0] + 
                  ['Generated 613'] * v_generated_data_613.shape[0] + ['Generated 352'] * v_generated_data_352.shape[0])

g_combined_data = np.vstack([existing_data, g_generated_data_612, g_generated_data_613, g_generated_data_352])
g_labels = np.array(['Existing'] * existing_data.shape[0] + ['Generated 612'] * g_generated_data_612.shape[0] + 
                  ['Generated 613'] * g_generated_data_613.shape[0] + ['Generated 352'] * g_generated_data_352.shape[0])

v_combined_data = np.vstack([existing_data, v_generated_data_flight, v_generated_data_ground])
v_labels = np.array(['Existing'] * existing_data.shape[0] + ['Generated Flight'] * v_generated_data_flight.shape[0] + 
                  ['Generated Ground'] * v_generated_data_ground.shape[0])

g_combined_data = np.vstack([existing_data, g_generated_data_flight, g_generated_data_ground])
g_labels = np.array(['Existing'] * existing_data.shape[0] + ['Generated Flight'] * v_generated_data_flight.shape[0] + 
                  ['Generated Ground'] * v_generated_data_ground.shape[0])

v_umap_model = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1)
v_umap_result = v_umap_model.fit_transform(v_combined_data)

g_umap_model = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1)
g_umap_result = v_umap_model.fit_transform(g_combined_data)

v_umap_df = pd.DataFrame(v_umap_result, columns=['UMAP Component 1', 'UMAP Component 2'])
v_umap_df['Data Type'] = (['Existing'] * existing_data.shape[0] + ['Generated 612'] * v_generated_data_612.shape[0] + 
                      ['Generated 613'] * v_generated_data_613.shape[0] + ['Generated 352'] * v_generated_data_352.shape[0])
g_umap_df = pd.DataFrame(g_umap_result, columns=['UMAP Component 1', 'UMAP Component 2'])
g_umap_df['Data Type'] = (['Existing'] * existing_data.shape[0] + ['Generated 612'] * g_generated_data_612.shape[0] + 
                  ['Generated 613'] * g_generated_data_613.shape[0] + ['Generated 352'] * g_generated_data_352.shape[0])

v_umap_df = pd.DataFrame(v_umap_result, columns=['UMAP Component 1', 'UMAP Component 2'])
v_umap_df['Data Type'] = (['Existing'] * existing_data.shape[0] + ['Generated Flight'] * v_generated_data_flight.shape[0] + 
                  ['Generated Ground'] * v_generated_data_ground.shape[0])
g_umap_df = pd.DataFrame(g_umap_result, columns=['UMAP Component 1', 'UMAP Component 2'])
g_umap_df['Data Type'] = (['Existing'] * existing_data.shape[0] + ['Generated Flight'] * g_generated_data_flight.shape[0] + 
                  ['Generated Ground'] * g_generated_data_ground.shape[0])

# colors = {'Existing': 'blue', 'Generated 612': 'yellow', 'Generated 613' : 'orange', 'Generated 352' : 'red'}
colors = {'Existing': 'blue', 'Generated Flight': 'orange', 'Generated Ground' : 'green'}

# Plot UMAP using Matplotlib
plt.figure(figsize=(8, 6))
for data_type in v_umap_df['Data Type'].unique():
    subset = v_umap_df[g_umap_df['Data Type'] == data_type]
    plt.scatter(subset['UMAP Component 1'], subset['UMAP Component 2'], label=data_type, color=colors[data_type], s=5, alpha=0.7)

plt.title('Vanilla CVAE UMAP: Color by Data Condition', fontsize=16)
plt.xlabel('UMAP Component 1', fontsize=14)
plt.ylabel('UMAP Component 2', fontsize=14)
plt.legend(title='Data Condition', fontsize=14)
plt.show()

plt.figure(figsize=(8, 6))
for data_type in g_umap_df['Data Type'].unique():
    subset = g_umap_df[g_umap_df['Data Type'] == data_type]
    plt.scatter(subset['UMAP Component 1'], subset['UMAP Component 2'], label=data_type, color=colors[data_type], s=10, alpha=0.7)

plt.title('GMA CVAE UMAP: Color by Data Condition', fontsize=16)
plt.xlabel('UMAP Component 1', fontsize=14)
plt.ylabel('UMAP Component 2', fontsize=14)
plt.legend(title='Data Condition', fontsize=14)
plt.show()


# DIFF GENE EXPRESSIOM

conditions = torch.tensor([0, 0, 14, 28, 1], dtype=torch.float32)
g_generated_data_flight_dif = van_model.generate(conditions, num_samples =10000).numpy()
conditions = torch.tensor([0, 0, 14, 28, 0], dtype=torch.float32)
g_generated_data_ground_dif = van_model.generate(conditions, num_samples = 10000).numpy()

var_names = adata.var_names
g_generated_data_flight_dif = AnnData(X = g_generated_data_flight_dif)
g_generated_data_flight_dif.obs['Flight'] = 1
g_generated_data_flight_dif.var_names = var_names
g_generated_data_ground_dif = AnnData(X = g_generated_data_ground_dif)
g_generated_data_ground_dif.obs['Flight'] = 0
g_generated_data_ground_dif.var_names = var_names

synth_flight_dif = concat([g_generated_data_flight_dif, g_generated_data_ground_dif], axis = 0)
synth_flight_dif.obs['Flight'] = synth_flight_dif.obs['Flight'].astype('category')
print(synth_flight_dif)

sc.tl.rank_genes_groups(synth_flight_dif, groupby='Flight', key_added='t-test')
sc.pl.rank_genes_groups(synth_flight_dif, n_genes=25, sharey=False, key="t-test")
#sc.pl.rank_genes_groups_heatmap(real_flight_dif_exp, n_genes=5, key='t-test', groupby='Flight', show_gene_labels=True)
sc.pl.rank_genes_groups_matrixplot(synth_flight_dif, n_genes=10, key="t-test", groupby="Flight")

real_flight_dif_exp = adata[adata.obs['Strain'] == 'B6129SF2/J']
real_flight_dif_exp.obs = real_flight_dif_exp.obs.drop(columns=['Strain', 'Sex', 'Age at Launch', 'Duration', 'dataset'])
real_flight_dif_exp.obs['Flight'] = real_flight_dif_exp.obs['Flight'].astype('category')

sc.tl.rank_genes_groups(real_flight_dif_exp, groupby='Flight', key_added='t-test')
sc.pl.rank_genes_groups(real_flight_dif_exp, n_genes=25, sharey=False, key="t-test")
#sc.pl.rank_genes_groups_heatmap(real_flight_dif_exp, n_genes=5, key='t-test', groupby='Flight', show_gene_labels=True)
sc.pl.rank_genes_groups_matrixplot(real_flight_dif_exp, n_genes=10, key="t-test", groupby="Flight")
print(real_flight_dif_exp)

synth_flight_genes = synth_flight_dif.uns['t-test']['names']
synth_flight_top_genes = [gene for sublist in synth_flight_genes[:25] for gene in sublist]

# Extract top differentially expressed genes for real flight
real_flight_genes = real_flight_dif_exp.uns['t-test']['names']
real_flight_top_genes = [gene for sublist in real_flight_genes[:25] for gene in sublist]
synth_genes_set = set(synth_flight_top_genes)
real_genes_set = set(real_flight_top_genes)

# Find common genes
common_genes = synth_genes_set.intersection(real_genes_set)

# Print the common genes
print("Common differentially expressed genes:", common_genes)

# MODEL STRUCTURES

# VANILLA STRUCT
print(van_model)

# GMA STRUCT
print(gma_model)