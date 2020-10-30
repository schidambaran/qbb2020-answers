#!/bin/bash

# Index mm10 genome using bismark
bismark_genome_preparation index/

# Move to folder containing methylation data
cd methylation_fastq/

# Map experiments to genome using bismark
bismark --genome ../index/ -1 SRR3083926_1.chr6.fastq -2 SRR3083926_2.chr6.fastq 
bismark --genome ../index/ -1 SRR3083929_1.chr6.fastq -2 SRR3083929_2.chr6.fastq

# Sort BAM files
samtools sort -o SRR3083926_sorted.bam SRR3083926_1.chr6_bismark_bt2_pe.bam
samtools sort -o SRR3083929_sorted.bam SRR3083929_1.chr6_bismark_bt2_pe.bam

# Index BAM files
samtools index SRR3083926_sorted.bam
samtools index SRR3083929_sorted.bam

# Extract methylation data
bismark_methylation_extractor --paired-end --comprehensive --bedGraph SRR3083926_1.chr6_bismark_bt2_pe.bam
bismark_methylation_extractor --paired-end --comprehensive --bedGraph SRR3083929_1.chr6_bismark_bt2_pe.bam
