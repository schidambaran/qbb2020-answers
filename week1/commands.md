## Question 1. Coverage Analysis

##### a. How long is the reference genome?

233,806 bp

```
samtools faidx ref.fa
cat ref.fa.fai
```

##### b. How many reads are provided and how long are they?

frag180.1.fq: 35,178 reads that are 100 bp long  
frag180.2.fq: 35,178 reads that are 100 bp long  
jump2k.1.fq: 70,355 reads that are 50 bp long  
jump2k.2.fq: 70,355 reads that are 50 bp long

```
fastqc frag180.1.fq
fastqc frag180.2.fq
fastqc jump2k.1.fq
fastqc jump2k.2.fq
```

##### c. How much coverage do you expect to have?

About 30x coverage for both frag180 and jump2k (both are paired-end reads)

frag180: ((35,178 * 100) / 233,806) * 2  
jump2k: ((70,355 * 50) / 233,806) * 2

##### d. Plot the average quality value across the length of the reads.

frag180.1.fq  
![frag180.1](https://github.com/schidambaran/qbb2020-answers/blob/master/week1/frag180.1.png)
frag180.2.fq  
![frag180.2](https://github.com/schidambaran/qbb2020-answers/blob/master/week1/frag180.2.png)
jump2k.1.fq  
![jump2k.1](https://github.com/schidambaran/qbb2020-answers/blob/master/week1/jump2k.1.png)
jump2k.2.fq  
![jump2k.2](https://github.com/schidambaran/qbb2020-answers/blob/master/week1/jump2k.2.png)

## Question 2. Kmer Analysis

##### a. How many kmers occur exactly 50 times?

##### b. What are the top 10 most frequently occurring kmers?

##### c. What is the estimated genome size based on the kmer frequencies?

##### d. How well does the GenomeScope genome size estimate compare to the reference genome?

## Question 3. De novo assembly

##### a. How many contigs were produced?

##### b. What is the total length of the contigs?

##### c. What is the size of your largest contig?

##### d. What is the contig N50 size?

## Question 4. Whole Genome Alignment

##### a. What is the average identify of your assembly compared to the reference?

##### b. What is the length of the longest alignment?

##### c. How many insertions and deletions are in the assembly?

## Question 5. Decoding the insertion

##### a. What is the position of the insertion in your assembly?

##### b. How long is the novel insertion?

##### c. What is the DNA sequence of the encoded message?

##### d. What is the secret message?
