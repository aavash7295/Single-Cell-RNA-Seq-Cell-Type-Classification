import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.io import mmread

# =====================================================
# Load Data
# =====================================================

barcodes = pd.read_csv(
    "count_matrix_barcodes.tsv",
    header=None
)[0].tolist()

genes = pd.read_csv(
    "count_matrix_genes.tsv",
    header=None
)[0].tolist()

metadata = pd.read_csv("metadata.csv")

matrix = mmread(
    "count_matrix_sparse.mtx"
).T.tocsr()

# =====================================================
# Create AnnData Object
# =====================================================

metadata.index = barcodes

adata = ad.AnnData(
    X=matrix,
    obs=metadata,
    var=pd.DataFrame(index=genes)
)

adata.obs_names = barcodes
adata.var_names_make_unique()

print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

# =====================================================
# Quality Control
# =====================================================

adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("Filtering low-quality cells...")

adata = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["total_counts"] >= 500) &
    (adata.obs["total_counts"] <= 25000) &
    (adata.obs["pct_counts_mt"] <= 15)
].copy()

print("Remaining cells:", adata.n_obs)

# =====================================================
# Normalization
# =====================================================

print("Normalizing data...")

sc.pp.normalize_total(
    adata,
    target_sum=10000
)

sc.pp.log1p(adata)

# =====================================================
# Highly Variable Genes
# =====================================================

print("Finding highly variable genes...")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True
)

print("Highly variable genes:", adata.n_vars)

# =====================================================
# PCA
# =====================================================

print("Running PCA...")

sc.pp.pca(
    adata,
    n_comps=50
)

# =====================================================
# Neighborhood Graph
# =====================================================

print("Computing neighbors...")

sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=50
)

# =====================================================
# UMAP
# =====================================================

print("Running UMAP...")

sc.tl.umap(adata)

# =====================================================
# Leiden Clustering
# Requires:
# pip install leidenalg igraph
# =====================================================

print("Running Leiden clustering...")

sc.tl.leiden(
    adata,
    resolution=0.5,
    flavor="igraph",
    directed=False,
    n_iterations=2
)

# =====================================================
# UMAP Cluster Plot
# =====================================================

sc.pl.umap(
    adata,
    color="leiden",
    save="_leiden_clusters.png",
    show=False
)

# =====================================================
# Marker Gene Discovery
# =====================================================

print("Finding marker genes...")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False,
    save="_marker_genes.png",
    show=False
)

# =====================================================
# Save Marker Genes
# =====================================================

marker_genes = pd.DataFrame(
    adata.uns["rank_genes_groups"]["names"]
)

marker_genes.head(20).to_csv(
    "marker_genes.csv",
    index=False
)

# =====================================================
# Temporary Cell Type Labels
# Replace after inspecting marker_genes.csv
# =====================================================

clusters = adata.obs["leiden"].unique()

cell_type_mapping = {
    cluster: f"Cell_Type_{cluster}"
    for cluster in clusters
}

adata.obs["cell_type"] = (
    adata.obs["leiden"]
    .map(cell_type_mapping)
)

# =====================================================
# Annotated UMAP
# =====================================================

sc.pl.umap(
    adata,
    color="cell_type",
    save="_cell_types.png",
    show=False
)

# =====================================================
# Save Processed Data
# =====================================================

adata.write(
    "processed_single_cell_data.h5ad"
)

print("\nAnalysis Complete!")
print("Files generated:")
print("- marker_genes.csv")
print("- processed_single_cell_data.h5ad")
print("- figures/umap_leiden_clusters.png")
print("- figures/rank_genes_groups_marker_genes.png")
print("- figures/umap_cell_types.png")
