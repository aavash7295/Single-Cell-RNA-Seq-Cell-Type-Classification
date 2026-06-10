# Single-Cell RNA-Seq Cell Type Classification

## Overview

This project implements a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline using Python and Scanpy. The workflow performs quality control, normalization, dimensionality reduction, clustering, marker gene identification, and cell type annotation on breast cancer single-cell transcriptomic data.

The goal of this project is to identify distinct cellular populations within the breast tumor microenvironment and classify them based on gene expression signatures.

---

## Objectives

* Process raw scRNA-seq count data
* Remove low-quality cells through quality control filtering
* Normalize and transform gene expression data
* Identify highly variable genes
* Perform dimensionality reduction using PCA and UMAP
* Cluster cells using the Leiden algorithm
* Discover cluster-specific marker genes
* Annotate cell populations using known biological markers
* Visualize cellular heterogeneity within the tumor microenvironment

---

## Technologies Used

* Python
* Scanpy
* AnnData
* Pandas
* NumPy
* SciPy
* Matplotlib

---

## Analysis Workflow

### 1. Data Loading

The pipeline imports:

* Gene expression count matrix
* Cell barcode information
* Gene annotations
* Cell metadata

### 2. Quality Control

Cells were filtered based on:

* Minimum genes detected: 200
* Minimum UMI counts: 500
* Maximum UMI counts: 25,000
* Maximum mitochondrial content: 15%

### 3. Normalization

Gene counts were normalized to 10,000 counts per cell and log-transformed.

### 4. Feature Selection

The top 2,000 highly variable genes were selected for downstream analysis.

### 5. Dimensionality Reduction

* Principal Component Analysis (PCA)
* Uniform Manifold Approximation and Projection (UMAP)

### 6. Clustering

Leiden clustering was performed to identify transcriptionally distinct cell populations.

### 7. Marker Gene Discovery

Differential expression analysis was performed using the Wilcoxon rank-sum test to identify cluster-specific marker genes.

### 8. Cell Type Annotation

Clusters were assigned biological identities based on known marker genes.

---

## Dataset Summary

| Metric                   | Value      |
| ------------------------ | ---------- |
| Initial Cells            | 100,064    |
| Cells After QC           | 89,150     |
| Highly Variable Genes    | 2,000      |
| Clustering Method        | Leiden     |
| Dimensionality Reduction | PCA + UMAP |

---

## Cell Type Annotations

| Cluster | Cell Type                          |
| ------- | ---------------------------------- |
| 0       | Endothelial Cells                  |
| 1       | Smooth Muscle / Pericytes          |
| 2       | B Cells / Plasma-like Cells        |
| 3       | Regulatory T Cells                 |
| 4       | Cytotoxic T / NK Cells             |
| 5       | NK Cells                           |
| 6       | Myeloid Cells                      |
| 7       | CD8 T Cells                        |
| 8       | Activated T Cells                  |
| 9       | T Cells                            |
| 10      | Fibroblasts                        |
| 11      | Plasma Cells                       |
| 12      | Luminal Epithelial / Tumor Cells   |
| 13      | Basal Epithelial Cells             |
| 14      | Luminal Epithelial / Tumor Cells   |
| 15      | Epithelial / Tumor Cells           |
| 16      | B Cells                            |
| 17      | Secretory Epithelial Cells         |
| 18      | Macrophages / Monocytes            |
| 19      | Inflammatory Macrophages           |
| 20      | Plasmacytoid Dendritic Cells       |
| 21      | Stress / Proliferating Cells       |
| 22      | Stromal Fibroblasts                |
| 23      | Secretory Epithelial Cells         |
| 24      | Luminal Epithelial / Tumor Cells   |
| 25      | Plasma Cells                       |
| 26      | Plasma Cells                       |
| 27      | Epithelial / Tumor Cells           |
| 28      | Basal-like Epithelial Cells        |
| 29      | Luminal Secretory Epithelial Cells |
| 30      | Basal Epithelial Cells             |
| 31      | Plasma Cells                       |

---

## Representative Marker Genes

### Endothelial Cells

* RAMP2
* PLVAP
* VWF
* CLDN5

### Fibroblasts

* DCN
* LUM
* COL1A1
* COL1A2

### B Cells

* CD19
* CD79A
* CD79B
* MS4A1

### Plasma Cells

* JCHAIN
* MZB1
* IGKC
* IGHG1

### Macrophages

* CD68
* LYZ
* FCGR3A
* CD14

### Basal Epithelial Cells

* KRT14
* KRT17
* KRT5
* KRT6A

### Luminal Tumor Cells

* AGR3
* ANKRD30A
* SLC39A6
* PIP

---

## Generated Outputs

```text
marker_genes.csv
figures/umap_leiden_clusters.png
figures/umap_cell_types.png
figures/rank_genes_groups_leiden_marker_genes.png
```

---

## Results

The analysis successfully identified diverse cellular populations within the breast cancer tumor microenvironment, including immune cells, stromal cells, endothelial cells, and multiple epithelial tumor cell populations. Marker gene analysis enabled biologically meaningful cell type annotation and visualization of tumor heterogeneity.

---

## Author

Aavash Pant

Bioinformatics | Computational Biology | Machine Learning | Genomics
