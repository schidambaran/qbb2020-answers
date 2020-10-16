#!/usr/bin/env python3

"""
Usage: ./week5_python.py
"""

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read in data on the number of CTCF binding sites present in functional regions of the mouse genome
binding_sites = pd.read_csv("CTCF_binding_sites.txt", sep = "\t")

# Read in data on the number of sites lost and gained during differentiation
binding_sites_lost = pd.read_csv("CTCF_binding_lost.bed", sep = "\t", header = None)
binding_sites_gain = pd.read_csv("CTCF_binding_gain.bed", sep = "\t", header = None)

fig, (ax1, ax2) = plt.subplots(ncols = 2)	

# Plot the number of CTCF binding sites in each region for each cell type
width = 0.4
x_vals = np.arange(3)

ax1.bar(x_vals, binding_sites["Number_of_Sites"][binding_sites["Cell"] == "G1E"], width = width, label = "G1E", color = "darkorange")
ax1.bar(x_vals + width, binding_sites["Number_of_Sites"][binding_sites["Cell"] == "ER4"], width = width, label = "ER4", color = "rebeccapurple")
ax1.set_xticks(x_vals + (width / 2))
ax1.set_xticklabels(["Exon", "Intron", "Promoter"])

ax1.set_title("CTCF Sites in Functional Regions")
ax1.set_xlabel("Region")
ax1.set_ylabel("Number of Sites")
ax1.legend(loc = "upper right")

# Plot the number of sites gained and lost during differentiation
ax2.bar(["Sites Lost", "Sites Gained"], [binding_sites_lost.shape[0], binding_sites_gain.shape[0]], color = ["steelblue", "maroon"])

ax2.set_title("CTCF Sites Lost or Gained\nAfter Differentiation")
ax2.set_ylabel("Number of Sites")

plt.tight_layout()
fig.savefig("CTCF_plots.png")
