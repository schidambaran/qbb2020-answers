#!/bin/bash

# Index the sacCer3 genome
bwa index sacCer3.fa 

# Align reads against sacCer3 genome
for SAMPLE in A01_09 A01_23 A01_27 A01_35 A01_62 A01_11 A01_24 A01_31 A01_39 A01_63
do
	bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -o ${SAMPLE}.sam sacCer3.fa ${SAMPLE}.fastq
done

# Create sorted BAM files for each sample
for SAMPLE in A01_09 A01_23 A01_27 A01_35 A01_62 A01_11 A01_24 A01_31 A01_39 A01_63
do
	samtools sort -O BAM -o ${SAMPLE}.bam ${SAMPLE}.sam 
done

# Identify genetic variants
freebayes -f sacCer3.fa --genotype-qualities -p 1 *.bam > sacCer3.vcf

# Filter variants based on genotype quality (only keep those whose estimated probability of being polymorphic is greater than 0.99)
vcffilter -f "QUAL > 20" sacCer3.vcf > sacCer3_sorted.vcf

# Decompose complex haplotypes
vcfallelicprimitives -k -g sacCer3_sorted.vcf > sacCer3_decomp.vcf

# Get yeast reference database
snpeff download R64-1-1.86

# Annotate VCF with predicted functional effects of variants
snpeff ann R64-1-1.86 sacCer3_decomp.vcf > sacCer3_ann.vcf