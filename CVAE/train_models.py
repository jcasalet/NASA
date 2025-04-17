import torch
import numpy as np
import src.utils as utils
import scanpy as sc
from scipy import sparse
import matplotlib.pyplot as plt
import pandas as pd
import src.vanilla_cvae as vanilla_cvae
import src.gma_cvae as gma_cvae

# To train the Vanilla CVAE, leave as 'vanilla', otherwise, change to gma
model_to_train = 'gma'

# Read in data
adata = sc.read_h5ad('data/corrected_data.h5ad')
adata.obs = adata.obs.drop(columns=['dataset'])
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]

# Create the data loaders
condition_cols = ['Strain', 'Sex', 'Age at Launch', 'Duration', 'Flight']
train_adata, val_adata = utils.split_adata(adata, condition_cols, val_split=0.3, random_state=42)
train_dataloader = utils.create_dataloader(train_adata, condition_cols, batch_size=128)
val_dataloader = utils.create_dataloader(val_adata, condition_cols, batch_size=128)

# GPU enable or disable
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Train/test vanilla cvae
if model_to_train == 'vanilla':
        # self.n_genes = n_genes
        # self.n_labels = n_labels
        # self.latent_size = latent_size
        # self.beta = kwargs.get('beta', 0.05)
        # self.dropout = kwargs.get('dropout', 0.03)
        # self.learn_rate = kwargs.get('lr', 0.001)
        # self.weight_decay = kwargs.get('wd', 0.001)
        # self.init_w = kwargs.get('init_w', False)
        # self.model = kwargs.get('model', 'trained_models/trained_vanilla_cvae.pt')
        # self.verbose = kwargs.get('verbose', True)
        # self.device = kwargs.get('device', torch.device('cpu'))
    van_cvae = vanilla_cvae.Vanilla_CVAE(n_genes=2000, n_labels=5, latent_size=64, beta=0.01, lr=0.001, wd=0.1, device=device)
    van_cvae.train_net(train_loader=train_dataloader, test_loader=val_dataloader, n_epochs=50)

# Train/test gma cvae
else:
    # Read in .gmt/create pathway mask
    gmt_dict = utils.read_gmt('data/BrainGMTv2_MouseOrthologs.gmt', min_g=0, max_g=2000)
    gm_mask = utils.create_pathway_mask(train_adata.var.index.tolist(), gmt_dict, n_labels=5, add_missing=5, fully_connected=True)
    print(gm_mask.shape)
        # self.pathway_mask = pathway_mask
        # self.n_pathways = self.pathway_mask.shape[1] - n_labels
        # self.n_genes = self.pathway_mask.shape[0]
        # self.n_labels = n_labels
        # self.beta = kwargs.get('beta', 0.05)
        # self.dropout = kwargs.get('dropout', 0.03)
        # self.learn_rate = kwargs.get('lr', 0.001)
        # self.weight_decay = kwargs.get('wd', 0.001)
        # self.init_w = kwargs.get('init_w', False)
        # self.model = kwargs.get('model', 'trained_models/trained_gma_cvae.pt')
        # self.verbose = kwargs.get('verbose', True)
        # self.device = kwargs.get('device', torch.device('cpu'))
    gma_cvae = gma_cvae.GMA_CVAE(n_labels=5, pathway_mask=gm_mask, beta=0.01, lr=0.001, wd=0.1, device=device)
    gma_cvae.train_net(train_loader = train_dataloader, test_loader = val_dataloader, n_epochs = 50)