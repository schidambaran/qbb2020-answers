#!/usr/bin/env python3

"""
Usage: ./week4_python.py
"""

# Import libraries
from fasta_parser import FASTAReader
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

# Read in protein alignments and separate out query alignment
protein_align = FASTAReader(open("aligned_seqs.fa"))
query_aa_align = protein_align[0]
protein_align = protein_align[1:]

# Read in query sequence
query_nt = FASTAReader(open("week4_query.fa"))[0]

# Read in BLAST output
blast_nt = FASTAReader(open("blast.fa"))

# Construt nucleotide alignment for query
query_nt_align = []
counter = 0

for aa in query_aa_align[1]:

	# If there is a gap in the protein alignment, add three gaps to the nucleotide alignment
	if (aa == "-"):
		query_nt_align.append("---")
	
	# Otherwise, add the corresponding codon to the nucleotide alignment
	else:
		query_nt_align.append(query_nt[1][counter:counter + 3])
		counter += 3

query_nt_align = "".join(query_nt_align)

# Construct nucleotide alignments for all other protein alignments
nt_align = []

for aa_seq, nt_seq in zip(protein_align, blast_nt):

	seq_align = []
	counter = 0

	for aa in aa_seq[1]:

		# If there is a gap in the protein alignment, add three gaps to the nucleotide alignment
		if (aa == "-"):
			seq_align.append("---")
		
		# Otherwise, add the corresponding codon to the nucleotide alignment
		else:
			seq_align.append(nt_seq[1][counter:counter + 3])
			counter += 3

	# Add nucleotide alignment to list
	nt_align.append("".join(seq_align))

# Create dictionary of codons and their corresponding amino acids
codontable = {"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M", 
"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", 
"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K", 
"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R", 
"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", 
"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", 
"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q", 
"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", 
"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", 
"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", 
"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E", 
"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", 
"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", 
"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L", 
"TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_", 
"TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W"}

# Go through all codons in query sequence and calculate their dN/dS ratios
codon_data = []
selection_diff = []

for base in range(0, len(query_nt_align), 3):

	codon = query_nt_align[base:base + 3]
	syn_change = 0
	non_change = 0
	ratio = 0.0

	# Skip "codons" that are just gaps
	if (codon == "---"):
		continue

	# Count the number of synonymous and nonsynonymous mutations in the codon in other alignments
	else:

		for seq in nt_align:

			new_codon = seq[base:base + 3]

			# Check if the new codon has changed and if both codons are in codontable
			if ((codon != new_codon) and (codon in codontable) and (new_codon in codontable)):

				# If both codons correspond to the same amino acid, add to synonymous mutations
				if (codontable[codon] == codontable[new_codon]):
					syn_change += 1

				# Otherwise, add to nonsynonymous changes
				else:
					non_change += 1

	if ((syn_change > 0) and (non_change > 0)):
		ratio = non_change / syn_change

	codon_data.append([codon, non_change, syn_change, ratio])

	# Add difference values (D) to list
	selection_diff.append(non_change - syn_change)

# Calculate standard deviation of all D values
sel_diff_std = scipy.stats.tstd(selection_diff)

# Calculate z score of each difference value
for entry in codon_data:

	D = entry[1] - entry[2]
	z_value = (0 - D) / sel_diff_std
	entry.append(z_value)

# Convert data into dataframe
codon_data = pd.DataFrame(codon_data)
codon_data.columns = ["Codon", "N", "S", "dN/dS", "z-score"]
codon_data = codon_data.loc[(codon_data["N"] != 0) & (codon_data["S"] != 0)]
codon_data["log2(dN/dS)"] = np.log2(codon_data["dN/dS"])

# Plot dN/dS ratios which have a p-value less than 0.05
# Significance is determined by comparing the z-score to -1.960
fig, ax = plt.subplots(figsize = (12, 8))

ax.scatter(codon_data.index, codon_data["log2(dN/dS)"])
ax.scatter(codon_data.index[codon_data["z-score"] < -1.960], codon_data["log2(dN/dS)"][codon_data["z-score"] < -1.960], color = "red")
plt.xticks(rotation = 45, ha = "right")

ax.set_title("dN/dS Ratios")
ax.set_xlabel("Codon Position")
ax.set_ylabel("log2(dN/dS)")

plt.tight_layout()
fig.savefig("week4_plot.png")	
