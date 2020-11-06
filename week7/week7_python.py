#!/usr/bin/env python3

"""
Usage: ./week7_python.py
"""

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
import seaborn as sns
import statsmodels.formula.api as smf

# Read in FPKM data
fpkm = pd.read_csv("all_annotated.csv")

# Filter dataset to genes with median expression greater than zero
fpkm = fpkm[fpkm.iloc[:, 2:].median(axis = 1) > 0]

# Transform data
fpkm = fpkm.iloc[:, :2].merge(np.log2(fpkm.iloc[:, 2:] + 0.1), left_index = True, right_index = True)

# Cluster data on genes
gene_clusters = linkage(fpkm.iloc[:, 2:])
gene_leaves = leaves_list(gene_clusters)

# Cluster data on samples
sample_clusters = linkage(fpkm.iloc[:, 2:].transpose())
sample_leaves = leaves_list(sample_clusters)

# Sort data by clusters
fpkm_clustered = fpkm.iloc[:, 2:].iloc[gene_leaves, sample_leaves]
fpkm_clustered.insert(loc = 0, column = "t_name", value = fpkm["t_name"].iloc[gene_leaves])
fpkm_clustered.insert(loc = 1, column = "gene_name", value = fpkm["gene_name"].iloc[gene_leaves])
fpkm_clustered = fpkm_clustered.set_index("t_name")

# Plot a heatmap of gene expression, clustered by both genes and samples
fig, ax = plt.subplots(figsize = (16, 12))

ax = sns.heatmap(data = fpkm_clustered.iloc[:, 2:])

ax.set_title("Drosophila Embryo Gene Expression")
ax.set_xlabel("Sample")
ax.set_ylabel("Transcript")

fig.savefig("heatmap.png")

# Create a dendrogram relating the samples to each other
fig, ax = plt.subplots(figsize = (16, 12))

ax = dendrogram(sample_clusters, labels = [list(fpkm_clustered.columns[1:])[i] for i in sample_leaves], leaf_rotation = 45)

plt.title("Dendrogram of Samples")

fig.savefig("dendrogram.png")
