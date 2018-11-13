
from crossmap.helpers import getBaseName
import os
import sys
from Bio import SeqIO

def extractTranscriptome(parsedArgs):
    for i in range(0,len(parsedArgs.annotations)):
        if parsedArgs.annotations[i].split(".")[-1] == "gtf":
            print("Annotation %s detected as gtf. Proceeding to transriptome extraction."%(os.path.basename(parsedArgs.annotations[i])))
            #get the transcriptome name
            transcriptome_name= getBaseName(parsedArgs.genomes[i])+"_transcriptome%s"%(i+1)+".fasta"

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
def readSimulation(parsedArgs, fasta_name,fasta_basename,file_number,read_len):
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


def simulateData(parsedArgs):
    file_num=0
    if parsedArgs.Simulation_type == "RNA":
        extractTranscriptome(parsedArgs)        
    for each_file in parsedArgs.fasta_names:
        fasta_basename = getBaseName(each_file)
        for rlen in parsedArgs.input_rlen:
            #print(each_file,fasta_basename,file_num,rlen)
            readSimulation(parsedArgs,each_file,fasta_basename,file_num,rlen)
            
        file_num+=1
    
    concateFastqFiles(parsedArgs)
        
def concateFastqFiles(parsedArgs):
    for rlen in parsedArgs.input_rlen:
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

