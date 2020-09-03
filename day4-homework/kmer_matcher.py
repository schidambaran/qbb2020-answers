#!/usr/bin/env python3

"""
Usage: ./kmer_matcher.py <target.fa> <query.fa> <k>
"""

# Load libraries
import sys
from fasta_iterator_class import FASTAReader

# Open files and get k-mer size
target = FASTAReader(open(sys.argv[1], "r"))
query = FASTAReader(open(sys.argv[2], "r"))
k = int(sys.argv[3])

kmers = {}

# Iterate through target file and get every k-mer
for tseq_id, tsequence in target:

	for i in range(0, len(tsequence) - k + 1):

		kmer = tsequence[i:i + k]
		kmers.setdefault(kmer, [])

		# Append sequence ID and k-mer position to dictionary (can have multiple tuples for one k-mer)
		kmers[kmer].append((tseq_id, i))

# Iterate through query sequence 
for qseq_id, qsequence in query:

	for j in range(0, len(qsequence) - k + 1):

		qkmer = qsequence[j:j + k]

		# If k-mer in query is also in target, print out information
		if qkmer in kmers:

			for entry in kmers[qkmer]:
				print("{}\t{}\t{}\t{}".format(entry[0], entry[1], j, qkmer))
