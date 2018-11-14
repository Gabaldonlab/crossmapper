import os
import math
from Bio import SeqIO
from crossmap.helpers import getBaseName





def prepareGenome(parsedArgs):
    genome_list=[]
    
    for i in range(0,len(parsedArgs.genomes)):
        genome_list.append(parsedArgs.genomes[i])    
    genome_concat = ' '.join(genome_list)
    
    
    cmd_genome_concat = f"cat {genome_concat} > {parsedArgs.out_dir}/concat.fasta"
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
    
    cmd_gtf_concat = f"cat {gtf_concat} > {parsedArgs.out_dir}/concat.gtf"
    print(cmd_gtf_concat)



###Mapping
    
def starIndex(parsedArgs):
    #calcualte concat.fasta genome size
    genome_len=0
    #for rec in SeqIO.parse(f"{parsedArgs.out_dir}/concat.fasta", 'fasta'):
     #   genome_len+=len(rec.seq)
        
    if genome_len > 3000000000:
        print("WARNING: concatenated genome size is larged than 3GB! " 
              +"\nMore than 30000 GB of RAM will be required for STAR mapping." )
    
    SA_index_size = min(14, round(math.log(genome_len,2)/2) - 1)
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


def bwaIndex(parsedArgs):
    algo = ""
        #calcualte concat.fasta genome size
    genome_len=0
    for rec in SeqIO.parse(f"{parsedArgs.out_dir}/concat.fasta", 'fasta'):
        genome_len+=len(rec.seq)
        
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
        
    cmd_bwa_index = "bwa index " \
    f"-p {parsedArgs.out_dir}/BWA_index/concat_BWA " \
    f"{algo} " \
    f"{parsedArgs.out_dir}/concat.fasta"

    print(cmd_bwa_index)
    print("Genome index for BWA is generated.")

    
    
def starMapping(parsedArgs,reads,rlen,read_layout):
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
f"--outFilterMultimapNmax 10000 " \
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
    
    
    
def bwaMapping(parsedArgs,reads,rlen,read_layout):
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
    

    
def mapping(parsedArgs):
    if parsedArgs.simulation_type == "RNA":
        starIndex(parsedArgs)
        for rlen in parsedArgs.input_rlen:
            se_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq")
            pe_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq") + " " + os.path.abspath(f"concat_{rlen}_read2.fastq")
            print(pe_mapping)
            if parsedArgs.read_layout == "SE":
                #cmd_remove_read2= "rm *_read2.fastq"
                starMapping(parsedArgs,se_mapping,rlen,parsedArgs.read_layout)
            elif parsedArgs.read_layout == "PE":
               starMapping(parsedArgs,pe_mapping,rlen,parsedArgs.read_layout)
            else:
                starMapping(parsedArgs,se_mapping,rlen,"SE")
                starMapping(parsedArgs,pe_mapping,rlen,"PE")
    else:
        bwaIndex(parsedArgs)
        for rlen in parsedArgs.input_rlen:
            se_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq")
            pe_mapping = os.path.abspath(f"concat_{rlen}_read1.fastq") + " " + os.path.abspath(f"concat_{rlen}_read2.fastq")
            if parsedArgs.read_layout == "SE":
               # cmd_remove_read2= "rm *_read2.fastq"
                bwaMapping(parsedArgs,se_mapping,rlen,parsedArgs.read_layout)
            elif parsedArgs.read_layout == "PE":
                bwaMapping(parsedArgs,pe_mapping,rlen,parsedArgs.read_layout)
            else:
                bwaMapping(parsedArgs,se_mapping,rlen,"SE")
                bwaMapping(parsedArgs,pe_mapping,rlen,"PE")