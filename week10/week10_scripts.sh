#!/bin/bash

# Reformat KRAKEN files for Krona graphs
for FILE in /Users/cmdb/qbb2020-answers/week10/week13_data/KRAKEN/*.kraken; do
	python week10_kraken_reformat.py ${FILE}
done

# Create Krona graphs
for KRAKEN_FILE in /Users/cmdb/qbb2020-answers/week10/week13_data/KRAKEN/*_edited.kraken; do

	# Get sample name
	SAMPLE=`echo ${KRAKEN_FILE} | cut -d '/' -f 8 | cut -d '_' -f 1`

	# Get filepath
	KRAKEN_PATH="/Users/cmdb/qbb2020-answers/week10/Krona/"${SAMPLE}

	# Create Krona graphs
	ktImportText -q -n root -o ${KRAKEN_PATH}.html ${KRAKEN_FILE}
done

# Index the assembly
bwa index week13_data/assembly.fasta

# Align reads against assembly
for SAMPLE in week13_data/READS/SRR492183 week13_data/READS/SRR492186 week13_data/READS/SRR492188 week13_data/READS/SRR492189 week13_data/READS/SRR492190 week13_data/READS/SRR492193 week13_data/READS/SRR492194 week13_data/READS/SRR492197; do
	bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -o ${SAMPLE}.sam week13_data/assembly.fasta ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
done

# Create sorted BAM files for each sample
for SAMPLE in week13_data/READS/SRR492183 week13_data/READS/SRR492186 week13_data/READS/SRR492188 week13_data/READS/SRR492189 week13_data/READS/SRR492190 week13_data/READS/SRR492193 week13_data/READS/SRR492194 week13_data/READS/SRR492197; do
	samtools sort -O BAM -o ${SAMPLE}.bam ${SAMPLE}.sam 
done

# Create depth file from BAM files
jgi_summarize_bam_contig_depths --outputDepth /Users/cmdb/qbb2020-answers/week10/metabat/depth.txt /Users/cmdb/qbb2020-answers/week10/week13_data/READS/*.bam

# Run MetaBAT to bin fragments
metabat2 -i /Users/cmdb/qbb2020-answers/week10/week13_data/assembly.fasta \
	-a /Users/cmdb/qbb2020-answers/week10/metabat/depth.txt \
	-o /Users/cmdb/qbb2020-answers/week10/metabat/bins

