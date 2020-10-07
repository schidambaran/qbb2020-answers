#!/bin/bash

# Combine BLAST output with query sequence
mv seqdump.txt blast.fa
cat week4_query.fa blast.fa > nt_seqs.fa

# Translate all sequences
transeq nt_seqs.fa aa_seqs.fa

# Align sequences
mafft aa_seqs.fa > aligned_seqs.fa