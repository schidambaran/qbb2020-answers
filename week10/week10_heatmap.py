#!/usr/bin/env python3

"""
Usage: ./week10_heatmap.py
"""

# Import libraries
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Read in abundance data
abundances = pd.read_csv("abundance_table.tab", sep = "\t", index_col = "Genomic bins")

# Add species information to data (bins are not in numerical order, so species data is manually added)
abundances.loc["bin.1", "Species"] = "Staphylococcus haemolyticus"
abundances.loc["bin.2", "Species"] = "Leuconostoc citreum"
abundances.loc["bin.3", "Species"] = "Staphylococcus lugdunensis"
abundances.loc["bin.4", "Species"] = "Enterococcus faecalis"
abundances.loc["bin.5", "Species"] = "Cutibacterium avidum"
abundances.loc["bin.6", "Species"] = "Staphylococcus epidermidis"
abundances.loc["bin.7", "Species"] = "Staphylococcus aureus"
abundances.loc["bin.8", "Species"] = "Anaerococcus prevotii"

abundances = abundances.set_index("Species")

# Plot heatmap of abundance values
heatmap = sns.clustermap(abundances, cmap = "viridis", figsize = (8, 6))
heatmap.savefig("heatmap.png")
