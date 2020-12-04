#!/bin/bash

# Create genomic partition file
hifive fends --length genome/mm9.len --binned 100000 genome_partition.hdf5

# Get counts of interaction reads
hifive hic-data --matrix data/WT_100kb/raw_\*.mat genome_partition.hdf5 interaction_counts.txt

# Create project file
hifive hic-project -f 25 -n 25 -j 100000 interaction_counts.txt project_file.txt

# Normalize data
hifive hic-normalize express -f 25 -w cis project_file.txt