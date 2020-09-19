#!/usr/bin/env python3

"""
Usage: python day4-lunch-1.py
"""

# Load libraries
import pandas as pd

# Open GTF
gtf = open("/Users/cmdb/data/genomes/BDGP6.Ensembl.81.gtf", "r")

# Initialize list for GTF data
gtf_3R = []

# Iterate through GTF to get information on the protein coding genes on chromosome 3R
for line in gtf:

	gtf_info = line.split("\t")

	if (gtf_info[0] == "3R") and (gtf_info[2] == "gene") and ("gene_biotype \"protein_coding\"" in gtf_info[-1]):
		
		# Add gene name, start position, and end position to list as a tuple
		gene_name = gtf_info[-1].split("; ")[2].split(" ")[1].strip("\"")
		gtf_3R.append((gene_name, int(gtf_info[3]), int(gtf_info[4])))

# Close GTF
gtf.close()

# Search for protein coding gene closest to position 21,378,950 on chromosome 3R
pos = 21378950

# Initialize variables for binary search
low = 0
high = len(gtf_3R) - 1
mid = 0
counter = 0
distance = 1000000000
gene_close = ""

while low <= high:

	counter += 1

	# Compare position to gene at mid
	mid = (low + high) // 2

	# If position is to the left of the gene at mid, move left
	if (pos < gtf_3R[mid][1]):
		high = mid - 1

		if (gtf_3R[mid][1] - pos) < distance:
			distance = gtf_3R[mid][1] - pos
			gene_close = gtf_3R[mid][0]

	# If position is to the right of the gene at mid, move right
	elif (pos > gtf_3R[mid][2]):
		low = mid + 1

		if (pos - gtf_3R[mid][2] < distance):
			distance = pos - gtf_3R[mid][2]
			gene_close = gtf_3R[mid][0]

	# Otherwise, position is in the gene at mid and distance is 0
	else:
		distance = 0
		gene_close = gtf_3R[mid][0]
		break

print("The closest protein coding gene to 3R:21,378,950 is {}, {} base pairs away. It took {} iterations to find this gene.".format(gene_close, distance, counter))
