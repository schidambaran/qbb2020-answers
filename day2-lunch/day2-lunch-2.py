#!/usr/bin/env python3

"""
Usage: ./day2-lunch-2.py fly_mapping.txt ~/data/results/stringtie/SRR072893/t_data.ctab
"""

# Import libraries
import sys

# Open files
mapping = open(sys.argv[1], "r")
c_tab = open(sys.argv[2], "r")
output = open("id_mapping.txt", "w")

# Use the same header for the output file as from the original c_tab file
output.write(c_tab.readline())

# Create a dictionary from the mapping file
fly_mapping = dict()

for line in mapping:
	split_line = line.split()
	fly_mapping[split_line[0]] = split_line[1]

# Iterate through c_tab file
for data in c_tab:

	c_tab_id = data.split()[8]

	# Check if FlyBase ID has a UniProt mapping
	if c_tab_id in fly_mapping:
		new_data = data.split()
		new_data[8] = fly_mapping[c_tab_id]
		new_data = "\t".join(new_data)
	# else:
	# 	new_data = data.split()
	# 	new_data[8]

# Close files
mapping.close()
c_tab.close()
output.close()
