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
from statsmodels.stats import multitest

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

# Reformat gene expression data for OLS and add columns for developmental stage and sex
fpkm_reformat = pd.melt(fpkm, id_vars = fpkm.columns[:2], value_vars = fpkm.columns[2:], value_name = "fpkm")
fpkm_reformat["stage"] = fpkm_reformat["variable"].str.split("_", expand = True)[1]
fpkm_reformat = fpkm_reformat.replace(["14A", "14B", "14C", "14D"], "14")
fpkm_reformat["stage"] = fpkm_reformat["stage"].astype(float)
fpkm_reformat["sex"] = fpkm_reformat["variable"].str.split("_", expand = True)[0]
fpkm_reformat.drop(columns = "variable", inplace = True)

# Use ordinary least squares (OLS) to test for genes that are differentially expressed across stages
p_values = []

for transcript in fpkm_reformat["t_name"].unique():

	# Get all expression data for one transcript
	transcript_data = fpkm_reformat[fpkm_reformat["t_name"] == transcript]

	# Use OLS to test if transcript is differentially expressed across stages
	model = smf.ols(formula = "fpkm ~ stage", data = transcript_data)
	results = model.fit()
	p_values.append([transcript, float(results.pvalues["stage"])])

# Convert p-value data to dataframe and calculate -log10(p-values)
p_values_df = pd.DataFrame(p_values, columns = ["Transcript", "p_values"])
p_values_df = p_values_df.sort_values(by = "p_values")
p_values_df["log_p_values"] = -1 * np.log10(p_values_df["p_values"])

# Calculate expected p-values
p_values_df["uniform_points"] = range(0, len(p_values_df))
p_values_df["uniform_pval"] = (p_values_df["uniform_points"] + 1) / len(p_values_df)
p_values_df["uniform_logP"] = -1 * np.log10(p_values_df["uniform_pval"])

# Generate QQ plot for p-values
fig, ax = plt.subplots()

ax.scatter(p_values_df["uniform_logP"], p_values_df["log_p_values"])
ax.plot([8, 0], [8, 0], color = "black")

ax.set_title("QQ Plot")
ax.set_xlabel("Expected -log10(p-value)")
ax.set_ylabel("Observed -log10(p-value)")

fig.savefig("qq_plot.png")

# Identify transcripts that are differential expressed at a 10% false discovery rate
p_values_df["fdr_0.10"] = multitest.multipletests(p_values_df["p_values"], method = "fdr_bh", alpha = 0.10)[0]

# Write these transcripts to an output file
p_values_df["Transcript"][p_values_df["fdr_0.10"]].to_csv("diff_expression.txt", index = False)

# Repeat analysis, but with sex as a covariate
p_values_cov = []
for transcript in fpkm_reformat["t_name"].unique():

	# Get all expression data for one transcript
	transcript_data = fpkm_reformat[fpkm_reformat["t_name"] == transcript]

	# Use OLS to test if transcript is differentially expressed across stages while controlling for sex
	model = smf.ols(formula = "fpkm ~ stage + sex", data = transcript_data)
	results = model.fit()
	p_values_cov.append([transcript, float(results.pvalues["stage"])])

# Convert p-value data to dataframe and calculate -log10(p-values)
p_values_cov_df = pd.DataFrame(p_values_cov, columns = ["Transcript", "p_values"])
p_values_cov_df = p_values_cov_df.sort_values(by = "p_values")
p_values_cov_df["log_p_values"] = -1 * np.log10(p_values_cov_df["p_values"])

# Calculate expected p-values
p_values_cov_df["uniform_points"] = range(0, len(p_values_cov_df))
p_values_cov_df["uniform_pval"] = (p_values_cov_df["uniform_points"] + 1) / len(p_values_cov_df)
p_values_cov_df["uniform_logP"] = -1 * np.log10(p_values_cov_df["uniform_pval"])

# Generate QQ plot for p-values
fig, ax = plt.subplots()

ax.scatter(p_values_cov_df["uniform_logP"], p_values_cov_df["log_p_values"])
ax.plot([8, 0], [8, 0], color = "black")

ax.set_title("QQ Plot with Sex as a Covariate")
ax.set_xlabel("Expected -log10(p-value)")
ax.set_ylabel("Observed -log10(p-value)")

fig.savefig("qq_plot_cov.png")

# Identify transcripts that are differential expressed at a 10% false discovery rate
p_values_cov_df["fdr_0.10"] = multitest.multipletests(p_values_cov_df["p_values"], method = "fdr_bh", alpha = 0.10)[0]

# Write these transcripts to an output file
p_values_cov_df["Transcript"][p_values_cov_df["fdr_0.10"]].to_csv("diff_expression_cov.txt", index = False)

# Calculate percent overlap in differentially expressed genes with and without sex as a covariate
transcripts_wo_cov = list(p_values_df["Transcript"][p_values_df["fdr_0.10"]])
transcripts_w_cov = list(p_values_cov_df["Transcript"][p_values_cov_df["fdr_0.10"]])
overlap = [t for t in transcripts_wo_cov if t in transcripts_w_cov]

percent_overlap = (len(overlap) / len(transcripts_wo_cov)) * 100

# Write percent overlap to output file
overlap_output = open("percent_overlap.txt", "w")
overlap_output.write("The percent overlap with and without sex as a covariate is " + str(percent_overlap) + "%.")
overlap_output.close()
