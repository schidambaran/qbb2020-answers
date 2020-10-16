#!/bin/bash

# Create chromosome 19 index
bowtie2-build bowtie_index/chr19.fa bowtie_index/chr19

# Map reads to reference genome (chromosome 19)
for sample in CTCF_ER4 CTCF_G1E input_ER4 input_G1E
do
  bowtie2 -x bowtie_index/chr19 -U g1e/${sample}.fastq -S mapped_reads/${sample}.sam -p 6
  samtools view -bSo mapped_reads/${sample}.bam mapped_reads/${sample}.sam
  samtools sort mapped_reads/${sample}.bam -o mapped_reads/${sample}.sorted.bam
  samtools index mapped_reads/${sample}.sorted.bam
done

# Determine chromosome size (62660198)
tail -n +2 bowtie_index/chr19.fa | wc -c
 
# Call ChIP-seq peaks
macs2 callpeak -t mapped_reads/CTCF_ER4.bam -c mapped_reads/input_ER4.bam --format=BAM --name=ER4 --gsize=63000000
macs2 callpeak -t mapped_reads/CTCF_G1E.bam -c mapped_reads/input_G1E.bam --format=BAM --name=G1E --gsize=63000000

# Identify sites where CTCF binding was lost or gained during differentiation
bedtools intersect -a G1E_peaks.narrowPeak -b ER4_peaks.narrowPeak -wo -v > CTCF_binding_lost.bed
bedtools intersect -a ER4_peaks.narrowPeak -b G1E_peaks.narrowPeak -wo -v > CTCF_binding_gain.bed

# Count the number of CTCF binding sites that overlap with a functional region in the mouse genome
bedtools intersect -a Mus_musculus.GRCm38.94_features.bed -b G1E_peaks.narrowPeak ER4_peaks.narrowPeak -C -names G1E ER4 > feature_overlap.bed

# Count the number of CTCF binding sites in each region in each cell type
echo -e "Cell\tRegion\tNumber_of_Sites" > CTCF_binding_sites.txt
for cell in G1E ER4
do
	for region in exon intron promoter
	do
		num=`grep ${cell} feature_overlap.bed | grep ${region} | cut -f 8 | paste -sd+ - | bc`
		echo -e "${cell}\t${region}\t${num}"
	done
done >> CTCF_binding_sites.txt
