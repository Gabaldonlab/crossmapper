
from crossmapper.helpers import getBaseName, getLogger
import os
import sys
from Bio import SeqIO
import subprocess
import crossmapper
import crossmapper.externalExec
import random




def renameChromosomes(parsedArgs):
    for i in range(0,len(parsedArgs.genomes)):
        with open(f"{parsedArgs.genomes[i]}", "r+") as original_fasta, open(f"{parsedArgs.chr_rename_fasta[i]}", "w") as renamed_fasta:
            for line in original_fasta.readlines():
                if line.startswith(">"):
                    line=line.rstrip()
                    line=line.replace(">", f">{parsedArgs.speciesPrefix[i]}_")
                    #print(line)
                    renamed_fasta.write(f"{line}\n")
                else:
                    line = line.rstrip()
                    renamed_fasta.write(f"{line}\n")
    
         
          
        if parsedArgs.simulation_type == "RNA": 
            with open(f"{parsedArgs.annotations[i]}","r+") as original_gff, open(f"{parsedArgs.chr_rename_gff[i]}", "w") as renamed_gff:
                for line in original_gff.readlines():
                    if not line.startswith("#"):
                        line=line.rstrip().split("\t")
                        line[0] = parsedArgs.speciesPrefix[i] + "_" + line[0]
                        line = "\t".join(line)
                        renamed_gff.write(f"{line}\n")
                    else:
                        renamed_gff.write(f"{line}\n")
            
    
    parsedArgs.genomes = parsedArgs.chr_rename_fasta
    if parsedArgs.simulation_type == "RNA":
        parsedArgs.annotations = parsedArgs.chr_rename_gff
    #print("NEW GENOMES",parsedArgs.genomes)            
     
def concatAnnotations(parsedArgs):
    logger = getLogger()
    if parsedArgs.simulation_type == "RNA" :
        ### concatenate gtf files
        gtf_list=[]
        for i in range(0,len(parsedArgs.genomes)):
            if parsedArgs.annotations[i].split(".")[-1] == "gtf":
                gtf_list.append(parsedArgs.annotations[i])
            else:
                gtf_name = getBaseName(parsedArgs.annotations[i]) + ".gtf"
                gtf_list.append(f"{parsedArgs.out_dir}/{gtf_name}")
        
        gtf_concat = ' '.join(gtf_list)
        
    #    cmd_gtf_concat = f"cat {gtf_concat} > {parsedArgs.out_dir}/concat.gtf"
        res = crossmapper.externalExec.execute(f"cat {gtf_concat}",
                                      "cat",
                                      f"{parsedArgs.out_dir}/concat.gtf",
                                      None,
                                      f"{parsedArgs.out_dir}")
        if not res.resCheck():
            sys.exit("Execution fail")
    parsedArgs.annotationsGTFConcat = f"{parsedArgs.out_dir}/concat.gtf"

def extractTranscriptome(parsedArgs):
    logger = getLogger()
    
    for i in range(0,len(parsedArgs.annotations)):
        if parsedArgs.annotations[i].split(".")[-1] == "gtf":
            logger.info("Annotation file %s detected as gtf. Proceeding to transriptome extraction."%(os.path.basename(parsedArgs.annotations[i])))
            #get the transcriptome name
            transcriptome_name = getBaseName(parsedArgs.genomes[i])+"_transcriptome%s"%(i+1)+".fasta"

            # extract the transcript
            
            cmd_gffread_extract = f"gffread " \
            f"-w {parsedArgs.out_dir}/{transcriptome_name} " \
            f"-g {parsedArgs.genomes[i]} " \
            f"{parsedArgs.annotations[i]}"
            
            #print(cmd_gffread_extract)
            res = crossmapper.externalExec.execute(cmd_gffread_extract,"gffreadExtract" , outDir = f"{parsedArgs.out_dir}")
            if not res.resCheck():
                sys.exit("Execution fail")
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] parsedArgs.annotations[i]
            logger.info("Transcriptome extracted for %s"%(os.path.basename(parsedArgs.genomes[i])))
            
        elif parsedArgs.annotations[i].split(".")[-1] == "gff":
			
            logger.info("Annotation file %s detected as gff. Converting to gtf using gffread."%(os.path.basename(parsedArgs.annotations[i])))
            
            #converting to gtf
            gtf_name = getBaseName(parsedArgs.annotations[i]) + ".gtf"
            
            cmd_gffread_convert = f"gffread " \
            f"{parsedArgs.annotations[i]} " \
            f"-T -o {parsedArgs.out_dir}/{gtf_name}"
            #print(cmd_gffread_convert)
            
            res = crossmapper.externalExec.execute(cmd_gffread_convert,"gffreadConvert" , outDir = f"{parsedArgs.out_dir}")
            if not res.resCheck(stdoutRemove=True,stdErrRemove = True):
                sys.exit("Execution fail")
            #gffread parsedArgs.annotations[i] -T -o gtf_name
            
            logger.info("GFF --> GTF conversion is done. Proceeding to transriptome extraction.")
            
            #get the transcriptome name
            transcriptome_name = getBaseName(parsedArgs.genomes[i])+"_transcriptome%s"%(i+1)+".fasta"
            
            
            cmd_gffread_extract = f"gffread " \
            f"-w {parsedArgs.out_dir}/{transcriptome_name} " \
            f"-g {parsedArgs.genomes[i]} " \
            f"{parsedArgs.out_dir}/{gtf_name}"
            
            #print(cmd_gffread_extract)
            res = crossmapper.externalExec.execute(cmd_gffread_extract,"gffreadExtract" , outDir = f"{parsedArgs.out_dir}")
            if not res.resCheck(stdoutRemove=True,stdErrRemove = True):
                sys.exit("Execution fail")
            # extract the transcript
            #gffread -w transcriptome_name -g parsedArgs.genomes[i] gtf_name
            logger.info("Transcriptome extracted for %s"%(os.path.basename(parsedArgs.genomes[i])))
        else:
            logger.error("Error: annotation file %s is neither gtf nor in gff. Please check the annotation file."%(os.path.basename(parsedArgs.annotations[i])))
            sys.exit("Execution Failed")



               
def readSimulation(parsedArgs, fasta_name,fasta_basename,file_number,read_len):
    logger = getLogger()
    fasta_len=0
    for rec in SeqIO.parse(f"{parsedArgs.genomeConcatFasta}", 'fasta'):
        fasta_len+=len(rec.seq)

    parsedArgs.simDir = os.path.join(parsedArgs.out_dir,"wgsim_output") 
    if os.path.isdir(f"{parsedArgs.simDir}") == False:
        logger.info(f"Creating {parsedArgs.simDir} directory.")
        os.makedirs(f"{parsedArgs.simDir}")
    ## if possible to assign, calculate N_reads, based on C, else use input value
    try:
        N_reads = round(parsedArgs.coverage[file_number]*fasta_len/read_len)
    except Exception as ex:
        N_reads = parsedArgs.N_read[file_number]
        
    random_seed=int(parsedArgs.random_seed) + random.randint(1,100000)
        
    cmd_wgsim = f"wgsim " \
f"-e {parsedArgs.error} " \
f"-d {parsedArgs.outer_dist} " \
f"-s {parsedArgs.s_dev} " \
f"-N {N_reads} " \
f"-1 {read_len} " \
f"-2 {read_len} " \
f"-r {parsedArgs.mut_rate} " \
f"-R {parsedArgs.indel_fraction} " \
f"-X {parsedArgs.indel_extend} " \
f"-S {random_seed} " \
f"-A {parsedArgs.discard_ambig} " \
f"{fasta_name} {parsedArgs.simDir}/{fasta_basename}_{read_len}_read1.fastq {parsedArgs.simDir}/{fasta_basename}_{read_len}_read2.fastq "
    
    crossmapper.externalExec.execute(cmd_wgsim,"cmd_wgsim", outDir = f"{parsedArgs.simDir}", overwrite = False)
    return cmd_wgsim


def simulateData(parsedArgs):
    
    
    
    if parsedArgs.simulation_type == "RNA":
        extractTranscriptome(parsedArgs)    
        concatAnnotations(parsedArgs)
    ## 
    for rlen in parsedArgs.input_rlen:
        parsedArgs.simulationOutputFiles[rlen] = []
    #for each_file in parsedArgs.fasta_names:
        #fasta_basename = getBaseName(each_file)
        #for rlen in parsedArgs.input_rlen:
        file_num=0
        for each_file in parsedArgs.fasta_names:
            
            fasta_basename = getBaseName(each_file)
            #print(each_file,fasta_basename,file_num,rlen)
            readSimulation(parsedArgs,each_file,fasta_basename,file_num,rlen)
            file_num+=1
        ## clean now
        concateFastqFiles(parsedArgs, rlen)
        
        
        
        
        
        
def concateFastqFiles(parsedArgs, rlen):
    logger = getLogger()
    
    #for rlen in parsedArgs.input_rlen:
    genome_list_r1=[]
    genome_list_r2=[]
    for i in range(0,len(parsedArgs.genomes)):
        if parsedArgs.simulation_type == "RNA":
            read_1 = parsedArgs.simDir + "/" + getBaseName(parsedArgs.genomes[i]) + "_transcriptome" + str(i+1) + "_" + str(rlen) + "_read1.fastq"
            genome_list_r1.append(read_1)
            
            read_2 = parsedArgs.simDir + "/" + getBaseName(parsedArgs.genomes[i]) + "_transcriptome" + str(i+1) + "_" + str(rlen) + "_read2.fastq"
            genome_list_r2.append(read_2)
            #print(genome_list_r2)
        else:
            read_1 = parsedArgs.simDir + "/" + getBaseName(parsedArgs.genomes[i]) + "_" + str(rlen) + "_read1.fastq"
            genome_list_r1.append(read_1)
                
            read_2 = parsedArgs.simDir + "/" + getBaseName(parsedArgs.genomes[i]) + "_" + str(rlen) + "_read2.fastq"
            genome_list_r2.append(read_2)
            
    genome_concat1=' '.join(genome_list_r1)
        #cmd_read1_concat = f"cat {genome_concat1} > {parsedArgs.out_dir}/concat_{rlen}_read1.fastq"
    res = crossmapper.externalExec.execute(f"cat {genome_concat1}",
                                      "cat", 
                                      f"{parsedArgs.simDir}/concat_{rlen}_read1.fastq",
                                      None,
                                      f"{parsedArgs.out_dir}")
    if not res.resCheck():
        sys.exit("Execution fail")

    
    genome_concat2=' '.join(genome_list_r2)
        
        # cmd_read2_concat = f"cat {genome_concat2} > {parsedArgs.simDir}/concat_{rlen}_read2.fastq"
    if parsedArgs.read_layout != "SE":
        res = crossmapper.externalExec.execute(f"cat {genome_concat2}",
                                          "cat",
                                          f"{parsedArgs.simDir}/concat_{rlen}_read2.fastq",
                                          None,
                                          f"{parsedArgs.out_dir}")
        if not res.resCheck():
            sys.exit("Execution fail")
#    else: ## no need ??
#        ## remove right reads files
#        try :
#            logger.debug(f"Removeing simulated reads 2 from wgsim {parsedArgs.out_dir}/concat_{rlen}_read2.fastq")
#            os.remove(f"{parsedArgs.out_dir}/concat_{rlen}_read2.fastq")
#        except:
#            logger.warning(f"Can not remove unwanted reads file  {parsedArgs.out_dir}/concat_{rlen}_read2.fastq")
    
    parsedArgs.simulationOutputFiles[rlen].append(f"{parsedArgs.simDir}/concat_{rlen}_read1.fastq")
    if parsedArgs.read_layout != "SE":
        parsedArgs.simulationOutputFiles[rlen].append(f"{parsedArgs.simDir}/concat_{rlen}_read2.fastq")

    ## cleanning temp files
    tmpFiles = []
    tmpFiles.extend(genome_list_r1)
    tmpFiles.extend(genome_list_r2)
    for tmpFile in tmpFiles:
        try:
            logger.debug(f"Deleteing tmp file {tmpFile}")
            os.remove(tmpFile)
            logger.debug(f"tmp file {tmpFile} delete")
        except Exception:
            logger.error(f"Can not delete tmp file {tmpFile}", exc_info=True)
            ## ignoe if before is ok
                
                
                
                