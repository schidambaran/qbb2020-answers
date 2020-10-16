#!/bin/bash

# Sort ER4 by lowest p-value (highest -log10(p-value)) to find 100 strongest peaks
sort -k8,8 -n -r ER4_peaks.narrowPeak | head -n 100 > ER4_100_peaks.narrowPeak

# Get sequences of peaks in FASTA format
bedtools getfasta -fi bowtie_index/chr19.fa -bed ER4_100_peaks.narrowPeak > ER4_100_peaks.fasta

# Find motifs in the 100 strongest peaks in the ER4 state
meme-chip -meme-maxw 20 -oc memechip -db motif_databases/JASPAR/JASPAR_CORE_2016.meme ER4_100_peaks.fasta

# Convert strongest motif from EPS to PDF
epstopdf memechip/meme_out/logo1.eps
cp memechip/meme_out/logo1.pdf .
