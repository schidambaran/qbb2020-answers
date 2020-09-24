#!/usr/bin/env python3

# Import libraries
import matplotlib.pyplot as plt

# Open VCF file
vcf = open("sacCer3_ann.vcf", "r")

# Initialize lists to hold data
read_depth = []
genotype_quality = []
allele_frequency = []
predicted_effects = {}

# Iterate through VCF file
for line in vcf:

	# Skip header
	if (line[0] == "#"):
		continue

	# Extract INFO and FORMAT data from each line
	line_split = line.split("\t")
	format_fields = line_split[9:]
	info_fields = line_split[7].split(";")

	# Iterate through information in FORMAT
	for format_data in format_fields:

		f_data = format_data.split(":")

		# Get read depth information from every sample
		if (f_data[2] != "."):
			read_depth.append(float(f_data[2]))

		# Get genotype quality information from every sample
		if (f_data[1] != "."):
			genotype_quality.append(float(f_data[1]))

	# Iterate through information in INFO
	for info_data in info_fields:

		# Get allele frequency of each variant
		if ("AF=" == info_data[0:3]):
			allele_freqs = info_data.replace("AF=", "").split(",")
			allele_freqs = [float(freq) for freq in allele_freqs]
			allele_frequency.extend(allele_freqs)

		# Get predicted effects of each variant
		if ("ANN=" == info_data[0:4]):
			annot = info_data.replace("ANN=", "").split(",")
			
			for data in annot:
				predicted_effects.setdefault(data.split("|")[1], 1)
				predicted_effects[data.split("|")[1]] += 1

# Close VCF file
vcf.close()

# Set up plots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize = (16, 10))

# Plot read depth distribution of variant genotypes
ax1.hist(read_depth, bins = 30, color = "teal")
ax1.set_title("Read Depths of Variant Genotypes")
ax1.set_xlabel("Read Depth")
ax1.set_ylabel("Frequency")
ax1.set_yscale("log")

# Plot quality distribution of variant genotypes
ax2.hist(genotype_quality, bins = 30, color = "goldenrod")
ax2.set_title("Quality of Variant Genotypes")
ax2.set_xlabel("Quality")
ax2.set_ylabel("Frequency")

# Plot allele frequency spectrum
ax3.hist(allele_frequency, bins = 30, color = "maroon")
ax3.set_title("Allele Frequencies of Variant Genotypes")
ax3.set_xlabel("Allele Frequency")
ax3.set_ylabel("Frequency")

# Plot summary of predicted effects of variants
ax4.bar(predicted_effects.keys(), predicted_effects.values(), color = "olivedrab")
ax4.set_title("Predicted Effects of Variant Genotypes")
ax4.set_xlabel("Predicted Effect")
ax4.set_ylabel("Frequency")
ax4.set_yscale("log")
plt.xticks(rotation = 45, ha = "right")

plt.tight_layout()

# Save plots
fig.savefig("week2_plots.png")
