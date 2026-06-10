import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import mmread

barcodes = pd.read_csv("count_matrix_barcodes.tsv", header=None)[0].tolist()
genes = pd.read_csv("count_matrix_genes.tsv", header=None)[0].tolist()
metadata = pd.read_csv("metadata.csv")
matrix = mmread("count_matrix_sparse.mtx").T.tocsr()

adata = ad.AnnData(X=matrix,obs=metadata, var=pd.DataFrame(index=genes))

#mitrocondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#filter
adata = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["total_counts"] >= 500) &
    (adata.obs["total_counts"] <= 25000) &
    (adata.obs["pct_counts_mt"] <= 15)
].copy()

#normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True)

#dimensionality reduction
sc.pp.pca(adata, n_comps = 50)
sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 50)
sc.tl.umap(adata)

#clustering
sc.tl.leiden(adata, resolution = 0.5)

#visualization
sc.pl.umap(adata, color = "leiden", save = "umap_leiden.png")

#gene marker discovery
sc.tl.rank_genes_groups(adata, groupby = "leiden", method = "wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes = 20, sharey = False, save = "rank_genes.png")

#maker gene table
marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
marker_genes.to_csv("marker_genes.csv", index=False)

#cell type annotation
cell_type_mapping = {
    '0': 'Cell Type A',
    '1': 'Cell Type B',
    '2': 'Cell Type C',
    # Add more mappings as needed
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_mapping)
sc.pl.umap(adata, color = "cell_type", save = "umap_cell_type.png")








