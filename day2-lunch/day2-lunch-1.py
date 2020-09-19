#!/usr/bin/env python3

"""
Usage: python day2-lunch-1.py
"""

# Open fly mapping file
fly = open("fly.txt", "r")

# Open output file
output = open("fly_mapping.txt", "w")
output.write("flybase_id\tuniprot_id\n")

# Iterate through file
for line in fly:
    
    # Check if line contains "DROME"
    if line.find("DROME") > -1:
        
        # Split line on whitespace
        split_line = line.rstrip().split()
        
        # Ignore lines with missing information
        if len(split_line) < 4:
            continue
            
        # Write FlyBase ID and Uniprot ID to output file
        output.write("{}\t{}\n".format(split_line[-1], split_line[-2]))

# Close files
fly.close()
output.close()