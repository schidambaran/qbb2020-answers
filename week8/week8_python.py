#!/usr/bin/env python3

"""
Usage: ./week8_python.py
"""

# Import libraries
import scanpy as sc

# Read 10x dataset
adata = sc.read_10x_h5("neuron_10k_v3_filtered_feature_bc_matrix.h5")

# Make variable names (in this case the genes) unique
adata.var_names_make_unique()

# Create a PCA plot before filtering the data
sc.tl.pca(adata)
sc.pl.pca(adata, title = "PCA Before Filtering", save = "_before_filtering.png")

# Create a PCA plot after filtering the data using the Zheng et al., 2017 approach
sc.pp.recipe_zheng17(adata, n_top_genes = 1000, log = True, plot = False, copy = False)
sc.tl.pca(adata)
sc.pl.pca(adata, title = "PCA After Filtering", save = "_after_filtering.png")

# Identify clusters in the data using the Leiden algorithm
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# Create a t-SNE plot showing the clusters
sc.tl.tsne(adata)
sc.pl.tsne(adata, color = "leiden", title = "Leiden Clusters t-SNE", save = "_leiden.png")

# Create a UMAP plot showing the clusters
sc.tl.umap(adata)
sc.pl.umap(adata, color = "leiden", title = "Leiden Clusters UMAP", save = "_leiden.png")

# Identify genes that distinguish each cluster using the t-test method
sc.tl.rank_genes_groups(adata, groupby = "leiden", method = "t-test")
sc.pl.rank_genes_groups(adata, groupby = "leiden", save = "_t_test.png")

# Identify genes that distinguish each cluster using the logistic regression method
sc.tl.rank_genes_groups(adata, groupby = "leiden", method = "logreg")
sc.pl.rank_genes_groups(adata, groupby = "leiden", save = "_logreg.png")
