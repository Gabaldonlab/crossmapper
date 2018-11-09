#!/usr/bin/python

'''

### Crossmap soft

'''

import os
import re
import sys
import argparse
import subprocess 



desc = """
  --
  crossmap.py Software
"""

soft_version = "0.1"


standard_rlen = [25, 50, 75, 100, 125, 150, 300]

###Create top level parser

#formatter_class=argparse.RawTextHelpFormatter - with this the text in help is free of formating
#formatter_class=argparse.ArgumentDefaultsHelpFormatter - this adds default values in help. Should be in all subparsers.

parser = argparse.ArgumentParser(prog = "crossmap.py",description = desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="sub_command help", dest = "Simulation_type")

parser_DNA = subparsers.add_parser("DNA", help = "Simulate DNA data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_DNA.add_argument("-specific to bwa", "--specific", type=int, default=5,
	help = "Help message")
	
parser_RNA = subparsers.add_parser("RNA", help = "Simulate RNA data", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_RNA.add_argument("-max_mismatch", "--outFilterMismatchNmax", type=int, default=10, metavar="Int",
	help = "From STAR manual: "
	+ " alignment will be output only if it has no more mismatches than this value")

parser.add_argument("-g", "--genomes",type=str, nargs=2, required=True,
	help="Specify the genome files in fasta format. Enter genome names separated by whitespace. "
	+ "\n NOTE: Keep the same order of listing for gtf/gff files")
	
parser.add_argument("-a", "--annotations",type=str, nargs=2, required=True,
	help="Specify the gtf/gff files. Enter the file names separated by whitespace. "
	+ "NOTE: Keep the same order of listing as for genome files")
	
parser.add_argument("-t", "--threads", type=int, default = 1,
	help = "Number of cores to be used for all multicore-supporting steps")
	

#### WGSIM parameters

parser.add_argument("-e", "--error", type=float, default=0.02,
	help = "Base error rate")
	
parser.add_argument("-d", "--outer_dist", type=int, default = 500,
	help = 	"Outer distance between the two reads. For example, in case of 2x50 reads, d=300 and s=0 "
			+ "the  mates will be 200 bp apart from each other.")

parser.add_argument("-s", "--s_dev", type=int, default = 30,
	help = "Standard deviation of outer distance.")
	
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-N", "--N_read", type = int, nargs=2,
    help = "The number of reads/read pairs to generate. This paremeter can not be used alongside with -C ")

group.add_argument("-C", "--coverage", type = int, nargs=2,
    help = "Generate the number of reads that reaches the specified coverage. Coverage is calculated as:"
		+ "C = N*rlen/L, " 
		+ "where L is the length of the genome/transcriptome")
	
parser.add_argument("-rlay", "--read_layout", type = str, choices=["SE","PE","both"], default = "SE",
    help = "Specify the read configuration - single-end (SE), paired-end (PE), or both (both)." 
		+ " If chosen 'both', the software will make separate analysis with each configuration")

parser.add_argument("-rlen", "--read_length", type=int, nargs = "*", default=50,
	help = "Specify the read length. Choose from the possible read lengths available for Illumina machines:"
		+ "25, 50, 75, 100, 125, 150, 300. The user can either enter a specific length, or specify a comma-separated"
		+ "list of desired read lengths. In the latter case, the software will perform the analysis for all specified"
		+ "values separatelly and will report mapping statistics in a form of a graph")
		
parser.add_argument("-r", "--mut_rate", type=float, default = 0.001,
	help = "Mutation rate.")

parser.add_argument("-R", "--indel_fraction", type=float, default = 0.015,
	help = "Fraction of indels.")
	
parser.add_argument("-X", "--indel_extend", type=float, default = 0.3,
	help = "Probability of an indel to be extended.")
	
parser.add_argument("-S", "--random_seed", type=int, default=(-1), 
	help = "Seed for random generator.")

parser.add_argument("-A", "--discard_ambig", type=float, default = 0.05,
	help = "Disgard if the fraction of ambiguous bases is higher than this number.")
	
parser.add_argument("-hapl", "--haplotype_mode", action = "store_true", default = False,
	help = "Haplotype mode. If specified, the haploid mutations will be simulated instead of diploid.")
	
parser.add_argument("-v", "--version", action = "version", \
    version = "%(prog)s \"v" + soft_version + "\"")


parsedArgs =  parser.parse_args()
#print(parsedArgs)
#~ print(len(parsedArgs.annotations))

####### Check if rlen numbers are correct
print(parsedArgs.read_length)

## check if there are duplicated lengths
if not len(set(parsedArgs.read_length)) == len(parsedArgs.read_length):
    sys.exit("Error: read lengths shoud not be duplicated!")    

## check if any length is not standard
for length in parsedArgs.read_length:
    #print(length)
    if not length in standard_rlen:
        sys.exit("Error: input read length %s is not a standard Illumina read length."%(length)
                 + "\nPlease refer to our help page (crossmap -h) to find standard read lengths.")
        
        
    


## If no arguments are given, just show the help and finish
if len(sys.argv) == 1:
   parser.print_help()
   sys.exit(1)

	
###TODO : check if genome files are in fasta format
    
    
def getBaseName(filename):
    if len(filename.split("."))>1:
        basename = '.'.join(filename.split(".")[0:-1])
    else:
        sys.exit("Error: please check the extensions of your input files."
                 +"\nGenome files should have .fa, .fasta or .fsa extensions."
                 +"\nGenome annotations should have .gtf or .gff extensions.")
    return basename
	
def extractTranscriptome():
    for i in range(0,len(parsedArgs.annotations)):
        if parsedArgs.annotations[i].split(".")[-1] == "gtf":
            print("Annotation %s detected as gtf. Continuing."%(parsedArgs.annotations[i]))
            #get the transcriptome name
            transcriptome_name=getBaseName(parsedArgs.genomes[i])+"_transcriptome_%s"%(i+1)+".fasta"

            # extract the transcript
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] parsedArgs.annotations[i]
            print("Transcriptome extracted for %s"%(parsedArgs.genomes[i]))
            
        elif parsedArgs.annotations[i].split(".")[-1] == "gff":
			
            print("Annotation file %s detected as gff. Converting to gtf using gffread."%(parsedArgs.annotations[i]))
            
            #converting to gtf
            gtf_name = getBaseName(parsedArgs.annotations[i])+".gtf"
            
            #gffread parsedArgs.annotations[i] -T -o gtf_name
            
            print("GFF --> GTF conversion is done")
            
            #get the transcriptome name
            transcriptome_name = getBaseName(parsedArgs.genomes[i])+"_transcriptome_%s"%(i+1)+".fasta"
            
            # extract the transcript
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] gtf_name
            print("Transcriptome extracted for %s"%(parsedArgs.genomes[i]))
        else:
            sys.exit("Error: annotation file %s is neither gtf nor in gff. Please check the annotation file."%(parsedArgs.annotations[i]))




fasta_names=[]
if parsedArgs.Simulation_type == "RNA":
    for i in range(0,len(parsedArgs.annotations)):
        transcriptome_name=getBaseName(parsedArgs.genomes[i])+"_transcriptome_%s"%(i+1)+".fasta"
        fasta_names.append(transcriptome_name)
else:
    for i in range(0,len(parsedArgs.genomes)):
        fasta_names.append(parsedArgs.genomes[i])
#print(fasta_names)


    


def readSimulation(fasta_list):
    if parsedArgs.Simulation_type == "RNA":
        extractTranscriptome()
        for rlen in input_rlen:
            print(rlen)
            if parsedArgs.read_layout == "SE":
                print("simulate data and remove mate2")
            elif parsedArgs.read_layout == "PE":
                print("simulate both mates")
            else:
                print("simulate both mates, map with both mates, and only with one mate")
            
# =============================================================================
#         for fasta in fasta_list():
#             cmd = """
#         wgsim \
# -e {0} \
# -d {1} \
# -s 30 \
# -N 10000000 \
# -1 50 \
# -2 50 \
# -r 0.001 \
# -R 0.01 \
# -X 0.1 \
# -S 134254 \
# -A 0.01 \
# input1.fasta org1_read1.fastq org1_read2.fastq
# """.format(parsedArgs.error)
#     print(cmd)
# 
# =============================================================================


readSimulation(fasta_names)

# =============================================================================
# if parsedArgs.Simulation_type == "RNA":
#     print("Simulating RNA")
#     extractTranscriptome()
# 	 
# else:
# 	print("Simulating DNA")
# 
# =============================================================================

#~ ###Simulate reads
#~ #input1

#~ wgsim \
#~ -e 0.02 \
#~ -d 300 \
#~ -s 30 \
#~ -N 10000000 \
#~ -1 50 \
#~ -2 50 \
#~ -r 0.001 \
#~ -R 0.01 \
#~ -X 0.1 \
#~ -S 134254 \
#~ -A 0.01 \
#~ input1.fasta org1_read1.fastq org1_read2.fastq


#~ #input2
#~ wgsim \
#~ -e 0.02 \
#~ -d 300 \
#~ -s 30 \
#~ -N 10000000 \
#~ -1 50 \
#~ -2 50 \
#~ -r 0.001 \
#~ -R 0.01 \
#~ -X 0.1 \
#~ -S 134254 \
#~ -A 0.01 \
#~ input1.fasta org2_read1.fastq org2_read2.fastq


#~ #Options: -e FLOAT      base error rate [0.020]
#~ #         -d INT        outer distance between the two ends [500] For example, with d=300 and s=0 in case of 2x100 reads, mates will be 100 bp apart from each other
#~ #         -s INT        standard deviation [50]
#~ #         -N INT        number of read pairs [1000000]
#~ #         -1 INT        length of the first read [70]
#~ #         -2 INT        length of the second read [70]
#~ #         -r FLOAT      rate of mutations [0.0010]
#~ #         -R FLOAT      fraction of indels [0.15]
#~ #         -X FLOAT      probability an indel is extended [0.30]
#~ #         -S INT        seed for random generator [-1]
#~ #         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#~ #         -h            haplotype mode



#~ ### Merge genomes and fastq files

#~ cat reference_genome1.fasta reference_genome2.fasta > concat_reference_genome.fasta

#~ cat org1_read1.fastq org2_read1.fastq > concat_read1.fastq
#~ cat org1_read2.fastq org2_read2.fastq > concat_read2.fastq

#~ #### Mapping for DNA with BWA-MEM

#~ #Indexing
#~ bwa index -p concat concat_reference_genome.fasta

#~ if gneomesize is larger than 3 GB , use -a bwtsw

#~ #Mapping
#~ bwa mem -t 6 concat concat_read1.fastq concat_read2.fastq | samtools sort -@6 -o concat.bam -

#~ samtools index concat.bam


#~ ### Mapping with STAR
#~ #Indexing
#~ mkdir STAR_index
#~ STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ./STAR_index --genomeFastaFiles concat_reference_genome.fasta --genomeSAindexNbases 10

#~ define  --genomeSAindexNbases as min(14, log2(GenomeLength)/2 - 1)


#~ #Mapping
#~ STAR --runThreadN 10 --genomeDir STAR_index --sjdbGTFfile GTF_file --sjdbOverhang 49 \
#~ --readFilesIn concat_read1.fastq concat_read2.fastq.gz \
#~ --readFilesCommand cat --outSAMtype BAM Unsorted --outFileNamePrefix  \
#~ --outTmpDir ~/TMP/TMPs --outFilterMismatchNmax 10 --outFilterMultimapNmax 10000

#~ samtools sort -@6 -o concat.bam Aligned.out.bam

#~ define --sjdbOverhang = rlen-1



######################################

#~ ### BWA PARAMETERS
#~ ## For indexing

#~ -p STR 	Prefix of the output database [same as db filename]
#~ -a STR 	Algorithm for constructing BWT index. Available options are:
#~ is 	IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.
#~ bwtsw 	Algorithm implemented in BWT-SW. This method works with the whole human genome.

#~ ## For mapping
#~ ### -t INT 	Number of threads [1] 
#~ -k INT 	Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]
#~ -w INT 	Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]
#~ -d INT 	Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [100]
#~ -r FLOAT 	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]
#~ -c INT 	Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]
#~ -P 	In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
#~ -A INT 	Matching score. [1]
#~ -B INT 	Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]
#~ -O INT 	Gap open penalty. [6]
#~ -E INT 	Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]
#~ -L INT 	Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [5]
#~ -U INT 	Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [9]
#~ ###-p 	Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details.
#~ -R STR 	Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. [null]
#~ -T INT 	Don’t output alignment with score lower than INT. This option only affects output. [30]
#~ ### this should be default -a 	Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
#~ #-C 	Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
#~ -H 	Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
#~ -M 	Mark shorter split hits as secondary (for Picard compatibility).
#~ -v INT 	Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. [3] 


#~ ### STAR parameters, For now I will use only basic parameters:

#~ ## For indexing

#~ ## For mapping









