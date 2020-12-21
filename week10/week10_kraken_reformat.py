#!/usr/bin/env python3

"""
Usage: ./week10_kraken_reformat.py
"""

# Import libraries
import sys

# Get filepath for output file
output_path = sys.argv[1].replace(".kraken", "_edited.kraken")

# Open new file for reformatted data
kraken_reformat = open(output_path, "w")

# Reformat KRAKEN file
with open(sys.argv[1], "r") as kraken_file:

	for line in kraken_file:

		# Remove first column
		edit_line = "".join(line.split("\t")[1:])

		# Replace semicolon delimiters with tabs
		tab_line = "\t".join(edit_line.split(";"))

		# Write reformatted line to output file
		kraken_reformat.write(tab_line)

# Close output file
kraken_reformat.close()