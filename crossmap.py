#!/usr/bin/python

'''

### Crossmap soft

'''

import os
import re
import sys
import argparse
import subprocess 
import math
from Bio import SeqIO



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
subparsers = parser.add_subparsers(help="Simulation type. Choose to simulate either DNA or RNA data", dest = "Simulation_type")
subparsers.required = True

parser_DNA = subparsers.add_parser("DNA",help = "Simulate DNA data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

parser.add_argument("-rlen", "--read_length", type=str, default="50",
	help = "Specify the read length. Choose from the possible read lengths available for Illumina machines:"
		+ "25,50,75,100,125,150,300. The user can either enter a specific length, or specify a COMMA-SEPARATED (!)(no spaces are allowed between commas)"
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
	

parser.add_argument("-o", "--out_dir", default = "crossmap_out", type = str,
                   help = "Specify the output directory for crossmap output files.")

parser.add_argument("-v", "--version", action = "version", \
    version = "%(prog)s \"v" + soft_version + "\"")


parsedArgs = parser.parse_args()

####### Check if rlen numbers are correct


###TODO: how to catch if there is a space????


## check if not all values can be converted to int
try:
    list(map(int,parsedArgs.read_length.split(",")))
except Exception:
    sys.exit("Hmm, seems you have strings or floats in read length values. Trying to find a bug?")

## convert list of strings to list of integers
input_rlen=list(map(int,parsedArgs.read_length.split(",")))
print(input_rlen)

## check if there are duplicated lengths
if not len(set(input_rlen)) == len(input_rlen):
    sys.exit("Error: read lengths shoud not be duplicated!")    

## check if any length is not standard
for length in input_rlen:
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
        
        basename = '.'.join(os.path.basename(filename).split(".")[0:-1])
    else:
        sys.exit("Error: please check the extensions of your input files."
                 +"\nGenome files should have .fa, .fasta or .fsa extensions."
                 +"\nGenome annotations should have .gtf or .gff extensions.")
    return basename


if os.path.isdir("./%s"%(parsedArgs.out_dir)) == True:
    print("%s already directory exists. Continuing."%(parsedArgs.out_dir))
else:
    cmd_mkdir = "mkdir ./%s"%(parsedArgs.out_dir)
    
    
fasta_names=[]
if parsedArgs.Simulation_type == "RNA":
    for i in range(0,len(parsedArgs.genomes)):
        transcriptome_name = getBaseName(parsedArgs.genomes[i]) + "_transcriptome%s"%(i+1) + ".fasta"
        fasta_names.append(os.path.abspath(transcriptome_name))
else:
    for i in range(0,len(parsedArgs.genomes)):
        fasta_names.append(os.path.abspath(parsedArgs.genomes[i]))
print(fasta_names)



def extractTranscriptome():
    for i in range(0,len(parsedArgs.annotations)):
        if parsedArgs.annotations[i].split(".")[-1] == "gtf":
            print("Annotation %s detected as gtf. Proceeding to transriptome extraction."%(os.path.basename(parsedArgs.annotations[i])))
            #get the transcriptome name
            transcriptome_name=getBaseName(parsedArgs.genomes[i])+"_transcriptome%s"%(i+1)+".fasta"

            # extract the transcript
            
            cmd_gffread_extract = f"gffread " \
            f"-w {parsedArgs.out_dir}/{transcriptome_name} " \
            f"-g {parsedArgs.genomes[i]} " \
            f"{parsedArgs.annotations[i]}"
            
            print(cmd_gffread_extract)
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] parsedArgs.annotations[i]
            print("Transcriptome extracted for %s"%(os.path.basename(parsedArgs.genomes[i])))
            
        elif parsedArgs.annotations[i].split(".")[-1] == "gff":
			
            print("Annotation file %s detected as gff. Converting to gtf using gffread."%(os.path.basename(parsedArgs.annotations[i])))
            
            #converting to gtf
            gtf_name = getBaseName(parsedArgs.annotations[i])+".gtf"
            
            cmd_gffread_convert = f"gffread " \
            f"{parsedArgs.annotations[i]} " \
            f"-T -o {parsedArgs.out_dir}/{gtf_name}"
            print(cmd_gffread_convert)
            
            #gffread parsedArgs.annotations[i] -T -o gtf_name
            
            print("GFF --> GTF conversion is done. Proceeding to transriptome extraction.")
            
            #get the transcriptome name
            transcriptome_name = getBaseName(parsedArgs.genomes[i])+"_transcriptome%s"%(i+1)+".fasta"
            
            
            cmd_gffread_extract = f"gffread " \
            f"-w {parsedArgs.out_dir}/{transcriptome_name} " \
            f"-g {parsedArgs.genomes[i]} " \
            f"{parsedArgs.out_dir}/{gtf_name}"
            
            print(cmd_gffread_extract)
            # extract the transcript
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] gtf_name
            print("Transcriptome extracted for %s"%(os.path.basename(parsedArgs.genomes[i])))
        else:
            sys.exit("Error: annotation file %s is neither gtf nor in gff. Please check the annotation file."%(os.path.basename(parsedArgs.annotations[i])))





#N_reads=59
def readSimulation(fasta_name,fasta_basename,file_number,read_len):
    fasta_len=0
    for rec in SeqIO.parse(f"{parsedArgs.out_dir}/concat.fasta", 'fasta'):
        fasta_len+=len(rec.seq)

    ## if possible to assign, calculate N_reads, based on C, else use input value
    try:
        N_reads = round(parsedArgs.coverage[file_number]*fasta_len/read_len)
    except Exception:
        N_reads = parsedArgs.N_read[file_number]
        
        
    wgsim_cmd = f"wgsim " \
f"-e {parsedArgs.error} " \
f"-d {parsedArgs.outer_dist} " \
f"-s {parsedArgs.s_dev} " \
f"-N {N_reads} " \
f"-1 {read_len} " \
f"-2 {read_len} " \
f"-r {parsedArgs.mut_rate} " \
f"-R {parsedArgs.indel_fraction} " \
f"-X {parsedArgs.indel_extend} " \
f"-S {parsedArgs.random_seed} " \
f"-A {parsedArgs.discard_ambig} " \
f"{fasta_name} {parsedArgs.out_dir}/{fasta_basename}_{read_len}_read1.fastq {parsedArgs.out_dir}/{fasta_basename}_{read_len}_read2.fastq"
    print(wgsim_cmd)
    return wgsim_cmd

        
        


def simulateData(fasta_list):
    file_num=0
    if parsedArgs.Simulation_type == "RNA":
        extractTranscriptome()        
    for each_file in fasta_list:
        fasta_basename = getBaseName(each_file)
        for rlen in input_rlen:
            #print(each_file,fasta_basename,file_num,rlen)
            readSimulation(each_file,fasta_basename,file_num,rlen)
            
        file_num+=1
    
                
simulateData(fasta_names)
                


## concatenate reference genomes
genome_list=[]

for i in range(0,len(parsedArgs.genomes)):
    genome_list.append(parsedArgs.genomes[i])    
genome_concat = ' '.join(genome_list)


cmd_genome_concat = f"{genome_concat} > {parsedArgs.out_dir}/concat.fasta"
print(cmd_genome_concat)

### concatenate gtf files
gtf_list=[]
for i in range(0,len(parsedArgs.genomes)):
    if parsedArgs.annotations[i].split(".")[-1] == "gtf":
        gtf_list.append(parsedArgs.annotations[i])
    else:
        gtf_name = getBaseName(parsedArgs.annotations[i]) + ".gtf"
        gtf_list.append(f"{parsedArgs.out_dir}/{gtf_name}")

gtf_concat = ' '.join(gtf_list)

cmd_gtf_concat = f"{gtf_concat} > {parsedArgs.out_dir}/concat.gtf"
print(cmd_gtf_concat)
    
    

## concatenate fastq files

for rlen in input_rlen:
    genome_list_r1=[]
    genome_list_r2=[]
    for i in range(0,len(parsedArgs.genomes)):
        if parsedArgs.Simulation_type == "RNA":
            read_1 = parsedArgs.out_dir + "/" + getBaseName(parsedArgs.genomes[i]) + "_transcriptome" + str(i+1) + "_" + str(rlen) + "_read1.fastq"
            genome_list_r1.append(read_1)
        
            read_2 = parsedArgs.out_dir + "/" + getBaseName(parsedArgs.genomes[i]) + "_transcriptome" + str(i+1) + "_" + str(rlen) + "_read2.fastq"
            genome_list_r2.append(read_2)
        #print(genome_list_r2)
        else:
            read_1 = parsedArgs.out_dir + "/" + getBaseName(parsedArgs.genomes[i]) + "_" + str(rlen) + "_read1.fastq"
            genome_list_r1.append(read_1)
            
            read_2 = parsedArgs.out_dir + "/" + getBaseName(parsedArgs.genomes[i]) + "_" + str(rlen) + "_read2.fastq"
            genome_list_r2.append(read_2)
        
    genome_concat1=' '.join(genome_list_r1)
    cmd_read1_concat = f"cat {genome_concat1} > {parsedArgs.out_dir}/concat_{rlen}_read1.fastq"
    print(cmd_read1_concat)
    
    
    genome_concat2=' '.join(genome_list_r2)
    cmd_read2_concat = f"cat {genome_concat2} > {parsedArgs.out_dir}/concat_{rlen}_read2.fastq"
    print(cmd_read2_concat)


###Mapping
    
def starIndex():
    #calcualte concat.fasta genome size
    genome_len=99999999999990
    #for rec in SeqIO.parse(f"{parsedArgs.out_dir}/concat.fasta", 'fasta'):
     #   genome_len+=len(rec.seq)
        
    if genome_len > 3000000000:
        print("WARNING: concatenated genome size is larged than 3GB! " 
              +"\nMore than 30 GB of RAM will be required for STAR mapping." )
    
    SA_index_size = min(14, round(math.log(26,2)/2) - 1)
    print("genomeSAindexNbases = %s"%(SA_index_size))
    print("Starting genome indexing with STAR.")
    if os.path.isdir(f"{parsedArgs.out_dir}/STAR_index") == True:
        print("STAR_index directory exists. Generating index files.")
    else:
        print(f"Creating {parsedArgs.out_dir}/STAR_index directory. Writing index files to STAR_index.")
        
        ###FOR SHELL
        #mkdir {parsedArgs.out_dir}/STAR_index
        
    cmd_star_index = "STAR " \
    f"--runThreadN {parsedArgs.threads} " \
    f"--genomeDir {parsedArgs.out_dir}/STAR_index " \
    f"--genomeFastaFiles {parsedArgs.out_dir}/concat.fasta " \
    f"--genomeSAindexNbases {SA_index_size}"
    print(cmd_star_index)
    print("Genome index for STAR is generated.")


def bwaIndex():
    algo = ""
        #calcualte concat.fasta genome size
    genome_len=0
   # for rec in SeqIO.parse(f"{parsedArgs.out_dir}/concat.fasta", 'fasta'):
       # genome_len+=len(rec.seq)
        
    if genome_len > 3000000000:
        print("Concatenated genome size is larged than 3GB. Using bwtsw algorithm for index generation" )
        algo = "-a bwtsw"
    
    print("Starting genome indexing with BWA.")
    if os.path.isdir(f"{parsedArgs.out_dir}/BWA_index") == True:
        print("BWA_index directory exists. Generating index files.")
    else:
        print(f"Creating {parsedArgs.out_dir}/BWA_index directory. Writing index files to BWA_index.")
        
        ### FOR SHELL
        #mkdir {parsedArgs.out_dir}/BWA_index
        #cd {parsedArgs.out_dir}/BWA_index
        
    cmd_bwa_index = "bwa index " \
    f"-p concat_BWA " \
    f"{algo} " \
    f"../concat.fasta"
    
    ### for shell
    #cd ..
    print(cmd_bwa_index)
    print("Genome index for BWA is generated.")

    
    
def starMapping(reads,rlen,read_layout):
    overhang = int(parsedArgs.read_length) - 1
    print("Starting STAR mapping.")
    cmd_star_mapping = "STAR " \
f"--runThreadN {parsedArgs.threads} " \
f"--genomeDir {parsedArgs.out_dir}/STAR_index " \
f"--sjdbGTFfile {parsedArgs.out_dir}/concat.gtf " \
f"--sjdbOverhang {overhang} " \
f"--readFilesIn {reads} " \
"--readFilesCommand cat --outSAMtype BAM Unsorted " \
f"--outFileNamePrefix concat_{rlen}_{read_layout}_ " \
f"--outFilterMismatchNmax {parsedArgs.outFilterMismatchNmax} " \
f"--outFilterMultimapNmax 10000" \
"--outTmpDir ~/TMP/TMPs"

    print(cmd_star_mapping)
    print("Mapping is finished. Started bam file sorting and indexing.")
    
    cmd_samtools_sort = "samtools sort " \
f"-@{parsedArgs.threads} " \
f"-o concat_{rlen}_{read_layout}_sorted.bam concat_{rlen}_{read_layout}_Aligned.out.bam"
    print(cmd_samtools_sort)
    print("Sorting is finished."
          +f"\nFinal bam file writen to concat_{rlen}_{read_layout}_sorted.bam")
    
    print("Starting bam indexing.")
    cmd_samtools_index = f"samtools index concat_{rlen}_{read_layout}_sorted.bam"
    print(cmd_samtools_index)
    print("Indexing is finished.")
    
    
    
def bwaMapping(reads,rlen,read_layout):
    print("Starting mapping with BWA.")
    cmd_bwa_mapping = "bwa mem " \
    f"-t {parsedArgs.threads} concat {reads} -a | " \
    f"samtools sort @{parsedArgs.threads} -o concat_{rlen}_{read_layout}_sorted.bam -"
    print(cmd_bwa_mapping)
    print("Mapping is finished."
          +f"\nFinal bam file writen to concat_{rlen}_{read_layout}_sorted.bam")
    
    print("Starting bam indexing.")
    cmd_samtools_index = f"samtools index concat_{rlen}_{read_layout}_sorted.bam"
    print(cmd_samtools_index)
    print("Indexing is finished.")
    

    
def mapping():
    if parsedArgs.Simulation_type == "RNA":
        starIndex()
        for rlen in input_rlen:
            se_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq")
            pe_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq") + " " + os.path.abspath(f"concat_{rlen}_read2.fastq")
            print(pe_mapping)
            if parsedArgs.read_layout == "SE":
                #cmd_remove_read2= "rm *_read2.fastq"
                starMapping(se_mapping,rlen,parsedArgs.read_layout)
            elif parsedArgs.read_layout == "PE":
               starMapping(pe_mapping,rlen,parsedArgs.read_layout)
            else:
                starMapping(se_mapping,rlen,"SE")
                starMapping(pe_mapping,rlen,"PE")
    else:
        bwaIndex()
        for rlen in input_rlen:
            se_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq")
            pe_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq") + " " + os.path.abspath(f"concat_{rlen}_read2.fastq")
            if parsedArgs.read_layout == "SE":
               # cmd_remove_read2= "rm *_read2.fastq"
                bwaMapping(se_mapping,rlen,parsedArgs.read_layout)
            elif parsedArgs.read_layout == "PE":
                bwaMapping(pe_mapping,rlen,parsedArgs.read_layout)
            else:
                bwaMapping(se_mapping,rlen,"SE")
                bwaMapping(pe_mapping,rlen,"PE")
            

mapping()








