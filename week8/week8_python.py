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

# Create UMAPs showing which clusters are enriched for certain genes and add info to dataset

sc.tl.leiden(adata, key_added = "clusters")

# Microglia (cluster 22)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Ccl4", save = "_Ccl4.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("22", "Microglia")

# Oligodendrocytes (cluster 17)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Olig2", save = "_Olig2.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("17", "Oligodendrocytes")

# Erythrocytes (cluster 14)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Hbb-bs", save = "_Hbb-bs.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("14", "Erythrocytes")

# Purkinje neurons (cluster 20)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Reln", save = "_Reln.png")

sc.tl.umap(adata)
sc.pl.umap(adata, color = "Lhx1", save = "_Lhx1.png")

adata.obs["clusters"] = adata.obs["clusters"].replace("20", "Purkinje Neurons")

# Endothelial cells (cluster 18)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Cldn5", save = "_Cldn5.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("18", "Endothelial Cells")

# Astrocytes (clusters 6 and 15)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Aldoc", save = "_Aldoc.png")

sc.tl.umap(adata)
sc.pl.umap(adata, color = "Ndrg2", save = "_Ndrg2.png")

adata.obs["clusters"] = adata.obs["clusters"].replace("6", "Astrocytes")
adata.obs["clusters"] = adata.obs["clusters"].replace("15", "Astrocytes")

# Intermediate progenitor cells (cluster 2)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "Eomes", save = "_Eomes.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("2", "Intermediate Progenitor Cells")

# Neural progenitor cells (cluster 10)
sc.tl.umap(adata)
sc.pl.umap(adata, color = "H2afz", save = "_H2afz.png")
adata.obs["clusters"] = adata.obs["clusters"].replace("10", "Neural Progenitor Cells")

sc.pl.umap(adata, color = "clusters", save = "_final.png")
