'''

### Crossmap soft

'''


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











