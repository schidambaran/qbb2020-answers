#!/usr/bin/env python3

"""
Usage: python day2-lunch-2.py fly_mapping.txt ~/data/results/stringtie/SRR072893/t_data.ctab ignore
"""

# Import libraries
import sys

# Open files
mapping = open(sys.argv[1], "r")
c_tab = open(sys.argv[2], "r")
output = open("id_mapping_ignore.txt", "w")

# Read in third argument
translation_arg = sys.argv[3]

# Use the same header for the output file as from the original c_tab file
output.write(c_tab.readline())

# Create a dictionary from the mapping file
fly_mapping = dict()

for line in mapping:
	split_line = line.split()
	fly_mapping[split_line[0]] = split_line[1]

# Limiting output to 100 lines
counter = 0

# Iterate through c_tab file
for data in c_tab:

	if (counter == 100):
		break

	c_tab_id = data.split()[8]

	# If FlyBase ID has a corresponding UniProt ID, replace the Flybase ID with the UniProt translation
	if c_tab_id in fly_mapping:
		new_data = data.split()
		new_data[8] = fly_mapping[c_tab_id]
		output.write("\t".join(new_data) + "\n")
		counter += 1

	# If no mapping exists, either replace the FlyBase ID with a default value or ignore that specific line
	else:

		# Ignore line
		if (translation_arg == "ignore"):
			continue

		# Replace FlyBase ID with .
		else:
			new_data = data.split()
			new_data[8] = "."
			output.write("\t".join(new_data) + "\n")
			counter += 1

# Close files
mapping.close()
c_tab.close()
output.close()
