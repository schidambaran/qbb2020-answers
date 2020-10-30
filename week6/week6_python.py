#!/usr/bin/env python3

"""
Usage: ./week6_python.py
"""

# Create dictionary for chromosome 6 genes and their coordinates
gene_coords = {}
with open("mm10_refseq_genes_chr6_50M_60M.bed", "r") as coord_file:
	for line in coord_file:
		data = line.split()
		gene_coords[data[12]] = (int(data[4]), int(data[5]))

# Open output file
output_file = open("methylation_enrichment.txt", "w")
output_file.write("Gene\tFold Change\n")


# Iterate through dictionary
for gene in gene_coords:

	E4_methyl = 0.0
	E5_5_methyl = 0.0

	# Calculate mean methylation signal in E4 cells
	with open("methylation_fastq/SRR3083926_1.chr6_bismark_bt2_pe.bedGraph", "r") as E4_file:

		# Skip header
		next(E4_file)

		# Find all methylation signals in this gene in E4 cells
		for data in E4_file:
			methyl_pos = int(data.split()[1])
			if (gene_coords[gene][0] <= methyl_pos <= gene_coords[gene][1]):
				E4_methyl += float(data.split()[3])

	# Calculate mean methylation signal of this gene in E4 cells
	E4_methyl /= (gene_coords[gene][1] - gene_coords[gene][0])

	# Calculate mean methylation signal in E5.5 cells
	with open("methylation_fastq/SRR3083929_1.chr6_bismark_bt2_pe.bedGraph", "r") as E5_5_file:

		# Skip header
		next(E5_5_file)

		# Find all methylation signals in this gene in E5.5 cells
		for data in E5_5_file:
			methyl_pos = int(data.split()[1])
			if (gene_coords[gene][0] <= methyl_pos <= gene_coords[gene][1]):
				E5_5_methyl += float(data.split()[3])

	# Calculate mean methylation signal of this gene in E5.5 cells
	E5_5_methyl /= (gene_coords[gene][1] - gene_coords[gene][0])

	# Calculate fold change in methylation signal from E4 to E5.5 cells and write to output file
	if (E4_methyl != 0):
		fold_change = (E5_5_methyl - E4_methyl) / E4_methyl
		output_file.write(str(gene) + "\t" + str(fold_change) + "\n")

# Close output file
output_file.close()
