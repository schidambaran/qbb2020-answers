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

`jellyfish count -m 21 -C -s 1000000 *.fq`

##### a. How many kmers occur exactly 50 times?

1,091 kmers occurred exactly 50 times

```
jellyfish histo mer_counts.jf > reads.histo
cat reads.histo
```

##### b. What are the top 10 most frequently occurring kmers?

GCCCACTAATTAGTGGGCGCC 105  
CCCACTAATTAGTGGGCGCCG 104  
CGCCCACTAATTAGTGGGCGC 104  
ACGGCGCCCACTAATTAGTGG 101  
AACAGGCCAGCTTATAAGCTG 98  
CAGGCCAGCTTATAAGCTGGC 98  
ACAGGCCAGCTTATAAGCTGG 97  
AGGCCAGCTTATAAGCTGGCC 95  
AGCATCGCCCACATGTGGGCG 83  
GCATCGCCCACATGTGGGCGA 82  

The number on the right is the frequency of each kmer.

`jellyfish dump -c mer_counts.jf | sort -k2,2nr | head -n 10`

##### c. What is the estimated genome size based on the kmer frequencies?

233,468 bp

##### d. How well does the GenomeScope genome size estimate compare to the reference genome?

The min "Genome Haploid Length" slightly underestimates the actual reference genome size; the estimate is 233,468 bp compared to the actual value of 233,806 bp.

## Question 3. De novo assembly

`spades.py --pe1-1 frag180.1.fq --pe1-2 frag180.2.fq --mp1-1 jump2k.1.fq --mp1-2 jump2k.2.fq -o asm -t 4 -k 31`

##### a. How many contigs were produced?

4

 `grep -c '>' asm/contigs.fasta`

##### b. What is the total length of the contigs?

234,467 bp

`cat asm/contigs.fasta.fai | cut -f 2 | paste -sd+ - | bc`

##### c. What is the size of your largest contig?

105,831 bp

`cat asm/contigs.fasta.fai | sort -k2,2nr`

##### d. What is the contig N50 size?

47,861 bp (calculated manually)

## Question 4. Whole Genome Alignment

##### a. What is the average identify of your assembly compared to the reference?

99.98%

```
dnadiff ref.fa asm/scaffolds.fasta 
cat out.report
```

##### b. What is the length of the longest alignment?

The longest alignment is 207,006 bp long in the reference and 206,998 bp long in the query.

```
nucmer ref.fa asm/scaffolds.fasta 
show-coords out.delta
```

##### c. How many insertions and deletions are in the assembly?

There is 1 insertion and 2 deletions in the assembly.

```
dnadiff ref.fa asm/scaffolds.fasta
cat out.report
```

## Question 5. Decoding the insertion

##### a. What is the position of the insertion in your assembly?

The insertion is between bases 206,998 and 207,711 in the assembly (on the reverse strand). In the reference, it's between bases 26,789 and 26,790.

`show-coords out.delta`

##### b. How long is the novel insertion?

The novel insertion is 711 bp long.

##### c. What is the DNA sequence of the encoded message?

>NODE_1_length_234497_cov_20.508040:206999-207710
TAACGATTTACATCGGGAAAGCTTAATGCAATTCACGCAGATATTCAGCTTAGAAGGTAC
GCAGCGGTGACGGGGTGCGGTCCATAATCTATGAAGCTATGAATTCGTACCTCAAGTAAT
GTTTTCTTCGCTGCAGTTCAGAAGTGATAAAGGTATCCCGCTTAGCCTGGCATACTTTGT
GCGTTCGTACCGCCCAGCATTAATGACTTGTGTAGGCAAGTAATGAACGACTCTTCTACG
CCGCGCCTAACCTCCGCACATAATGGCAGCATGTGGTAGTTACATACGCACAGAAGTGGT
TCGGTTTTAACTATAGTCAGATATGAATAAGCTGCGTGTGTCGTTGTGTCGGCGTGTCGT
ACTTACCTCCTGACATAGGTGAATTTCAGCCTACTGTAAGTTTGGAGTCGCGCTCTTTTC
TTATTATATTCTTTGGTATGTGTGTGATGGGTTCGGGCGTGTATTGATGTCTCTAAGGCT
CATGTTAGTGTTTATTTGGTCAGTTATGACGGTGTTCCTGTCGTACGTGTTGGCTTAGCG
GACTTGTAGACGGGATCAAGGTTGTCTGACCCTCCGGTCGACCGTGGGTCGGCCGTCCCG
GCCAGAATACAAGCCGCTTAGACTTTCGAAAGAGGGTAAGTTACTACGCGCGAACGTTAT
ACCTCGTTTCAGTATGCACTCCCTTAAGTCACTCAGAAAAGACTAAGGGGCT

```
samtools faidx asm/scaffolds.fasta NODE_1_length_234497_cov_20.508040:206999-207710 > encoded_message.fasta
cat encoded_message.fasta
```

##### d. What is the secret message?

The decoded message :  Congratulations to the 2020 CMDB @ JHU class!  Keep on looking for little green aliens...

`python ported_decoder.py --decode --rev_comp --input encoded_message.fasta`
