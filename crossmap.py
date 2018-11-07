'''

### Crossmap soft

'''


import os
import re
import sys
import argparse
import subprocess 


## Default values for certain variables
default_cores = 1
standard_rlen=["25", "50", "75", "100", "125", "150", "300"]

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--mode", type = str, choices=["DNA","RNA"], default = "DNA",
    help = "Specify whether to simulate DNA or RNA reads")

parser.add_argument("-G", "--genomes",type=str, nargs=2, required=True,
	help="Specify the genome files in fasta format. Enter genome names separated by whitespace. "
	+ "NOTE: Keep the same order of listing for gtf/gff files")
	
parser.add_argument("-A", "--annotations",type=str, nargs=2, required=True,
	help="Specify the gtf/gff files. Enter the file names separated by whitespace. "
	+ "NOTE: Keep the same order of listing for genome files")

parser.add_argument("-rlay", "--read_layout", type = str, choices=["SE","PE","both"], default = "SE"
    help = "Specify the read configuration - single-end (SE), paired-end (PE), or both (both)." 
		+ " If chosen 'both', the software will make separate analysis with each configuration")

parser.add_argument("-rlen", "--read_length", type=str, default="50",
	help = "Specify the read length. Choose from the possible read lengths available for Illumina machines:"
		+ "25, 50, 75, 100, 125, 150, 300. The user can either enter a specific length, or specify a comma-separated"
		+ "list of desired read lengths. In the latter case, the software will perform the analysis for all specified"
		+ "values separatelly and will report mapping statistics in a form of a graph")
		### we can check if input if correct with
		###	if set(rlen).issubset(standard_rlen)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-N", "--N_read", type = int, nargs=2,
    help = "The number of reads/read pairs to generate. This paremeter can not be used alongside with -C ")

group.add_argument("-C", "--coverage", type = int, nargs=2,
    help = "Generate the number of reads that reaches the specified coverage. Coverage is calculated as:"
		+ "C = N*rlen/L, " 
		+ "where L is the length of the genome/transcriptome")
		
parser.add_argument("-t", "--threads", type=int, default = 1,
	help = "Number of cores to be used for all multicore-supporting steps")
	





### Gffread command
#genome1
gffread -w output_transcripts1.fasta -g reference_genome1.fasta annotations1.gtf

#genome2

gffread -w output_transcripts2.fasta -g reference_genome2.fasta annotations2.gtf



###Simulate reads
#input1

wgsim \
-e 0.02 \
-d 300 \
-s 30 \
-N 10000000 \
-1 50 \
-2 50 \
-r 0.001 \
-R 0.01 \
-X 0.1 \
-S 134254 \
-A 0.01 \
input1.fasta org1_read1.fastq org1_read2.fastq


#input2
wgsim \
-e 0.02 \
-d 300 \
-s 30 \
-N 10000000 \
-1 50 \
-2 50 \
-r 0.001 \
-R 0.01 \
-X 0.1 \
-S 134254 \
-A 0.01 \
input1.fasta org2_read1.fastq org2_read2.fastq


#Options: -e FLOAT      base error rate [0.020]
#         -d INT        outer distance between the two ends [500] For example, with d=300 and s=0 in case of 2x100 reads, mates will be 100 bp apart from each other
#         -s INT        standard deviation [50]
#         -N INT        number of read pairs [1000000]
#         -1 INT        length of the first read [70]
#         -2 INT        length of the second read [70]
#         -r FLOAT      rate of mutations [0.0010]
#         -R FLOAT      fraction of indels [0.15]
#         -X FLOAT      probability an indel is extended [0.30]
#         -S INT        seed for random generator [-1]
#         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#         -h            haplotype mode



### Merge genomes and fastq files

cat reference_genome1.fasta reference_genome2.fasta > concat_reference_genome.fasta

cat org1_read1.fastq org2_read1.fastq > concat_read1.fastq
cat org1_read2.fastq org2_read2.fastq > concat_read2.fastq

#### Mapping for DNA with BWA-MEM

#Indexing
bwa index -p concat concat_reference_genome.fasta


#Mapping
bwa mem -t 6 concat concat_read1.fastq concat_read2.fastq | samtools sort -@6 -o concat.bam -

samtools index concat.bam


### Mapping with STAR
#Indexing
mkdir STAR_index
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ./STAR_index --genomeFastaFiles concat_reference_genome.fasta --genomeSAindexNbases 10

#Mapping
STAR --runThreadN 10 --genomeDir STAR_index --sjdbGTFfile GTF_file --sjdbOverhang 49 \
--readFilesIn concat_read1.fastq concat_read2.fastq.gz \
--readFilesCommand cat --outSAMtype BAM Unsorted --outFileNamePrefix  \
--outTmpDir ~/TMP/TMPs --outFilterMismatchNmax 10 

samtools sort -@6 -o concat.bam Aligned.out.bam











