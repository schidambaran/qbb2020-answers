#!/usr/bin/env python3

"""
Usage: ./week10_bin_classification.py
"""

# Create list of filepaths to binning results
files = ["bins/bin.1.fa", "bins/bin.2.fa", "bins/bin.3.fa", "bins/bin.4.fa", 
"bins/bin.5.fa", "bins/bin.6.fa", "bins/bin.7.fa", "bins/bin.8.fa"]

# Initialize output file
output = open("bin_classification.txt", "w")

# Iterate through all files
for filepath in files:

	# Open file
	with open(filepath, "r") as file:

		# Get name of each FASTA sequence
		fasta = []
		for line in file:
			if line[0] == ">":
				fasta.append(line.strip().lstrip(">"))

	# Read in KRAKEN scaffoled classification
	kraken_file = open("week13_data/KRAKEN/assembly.kraken", "r")
	kraken = kraken_file.readlines()
	kraken_file.close()

	# Initialize dictionaries to contain information on bin classification
	order = {}
	family = {}
	species = {}

	# Iterate through KRAKEN classification for matches to the bin results
	for kraken_data in kraken:

		if any(sample in kraken_data for sample in fasta):

			sample_name = kraken_data.split("\t")[0]
			sample_classification = kraken_data.split("\t")[1].split(";")

			# Extract KRAKEN classification data
			# Only look at lines with complete information
			if (len(sample_classification) > 9):

				species.setdefault(sample_classification[9], 0)
				species[sample_classification[9]] += 1

				family.setdefault(sample_classification[7], 0)
				family[sample_classification[7]] += 1

				order.setdefault(sample_classification[6], 0)
				order[sample_classification[6]] += 1


	output.write(filepath.rstrip(".fa").replace("bins/", "") + "\n\n")

	output.write("Species\tCount\n")
	for s_key, s_value in species.items():
		output.write(s_key + "\t" + str(s_value) + "\n")
	output.write("\n")

	output.write("Family\tCount\n")
	for f_key, f_value in family.items():
		output.write(f_key + "\t" + str(f_value) + "\n")
	output.write("\n")

	output.write("Order\tCount\n")
	for o_key, o_value in order.items():
		output.write(o_key + "\t" + str(o_value) + "\n")
	output.write("\n")

# Close output file
output.close()
