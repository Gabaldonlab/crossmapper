import os
import sys
import math
from Bio import SeqIO
from crossmapper.helpers import getBaseName, getLogger
import crossmapper



def concatGeneomes(parsedArgs):
    logger = getLogger()
    genome_list=[]
    for i in range(0,len(parsedArgs.chr_rename_fasta)):
        genome_list.append(parsedArgs.chr_rename_fasta[i])    
    genome_concat = ' '.join(genome_list)
    
    parsedArgs.genomeConcatFasta = f"{parsedArgs.out_dir}/concat.fasta"

#    cmd_genome_concat = f"cat {genome_concat} > {parsedArgs.out_dir}/concat.fasta"
    res = crossmapper.externalExec.execute(f"cat {genome_concat}",
                                  "cat",
                                  f"{parsedArgs.genomeConcatFasta}",
                                  None,
                                  f"{parsedArgs.out_dir}")
    if not res.resCheck():
        sys.exit("Execution fail")


###Mapping
    
def starIndex(parsedArgs):
    logger = getLogger()
    #calcualte concat.fasta genome size
    genome_len=0
    for rec in SeqIO.parse(parsedArgs.genomeConcatFasta, 'fasta'):
        genome_len+=len(rec.seq)
        
    if genome_len > 3000000000:
        logger.warning("Concatenated genome size is larged than 3 Gb! More than 30 GB of RAM will be required for STAR mapping." )
    
    SA_index_size = min(14, round(math.log(genome_len,2)/2) - 1)
    logger.debug("genomeSAindexNbases = %s"%(SA_index_size))
    logger.info("Starting genome indexing with STAR.")
    
    star_index = f"{parsedArgs.out_dir}/STAR_index"
    
    if os.path.isdir(f"{star_index}") == True:
        logger.info("STAR_index directory exists. Generating index files.")
    else:
        logger.info(f"Creating {star_index} directory. Writing index files to STAR_index.")
        os.makedirs(f"{star_index}")
        
        
        
    cmd_star_index = "STAR --runMode genomeGenerate " \
    f"--runThreadN {parsedArgs.threads} " \
    f"--genomeDir {star_index} " \
    f"--genomeFastaFiles {parsedArgs.genomeConcatFasta} " \
    f"--genomeSAindexNbases {SA_index_size}"
    
    res = crossmapper.externalExec.execute(cmd_star_index,"STAR_index" , outDir = f"{parsedArgs.out_dir}" )
    if not res.resCheck( stdoutRemove = True ):
        sys.exit("Execution Failed")
    logger.info("Genome index for STAR is generated.")
    parsedArgs.starIndex = star_index


def bwaIndex(parsedArgs):
    logger = getLogger()
    algo = ""
        #calcualte concat.fasta genome size
    genome_len=0
    for rec in SeqIO.parse(parsedArgs.genomeConcatFasta, 'fasta'):
        genome_len+=len(rec.seq)
        
    if genome_len > 3000000000:
        logger.info("Concatenated genome size is larged than 3GB. Using bwtsw algorithm for index generation" )
        algo = "-a bwtsw"
    
    logger.info("Starting genome indexing with BWA.")
    if os.path.isdir(f"{parsedArgs.out_dir}/BWA_index") == True:
        logger.info("BWA_index directory exists. Generating index files.")
    else:
        logger.info(f"Creating {parsedArgs.out_dir}/BWA_index directory. Writing index files to BWA_index.")
        os.makedirs(f"{parsedArgs.out_dir}/BWA_index")
        
    cmd_bwa_index = "bwa index " \
    f"-p {parsedArgs.out_dir}/BWA_index/concat_BWA " \
    f"{algo} " \
    f"{parsedArgs.genomeConcatFasta}"

    # print(cmd_bwa_index)
    res = crossmapper.externalExec.execute(cmd_bwa_index,"BWA_index", outDir = f"{parsedArgs.out_dir}")
    if not res.resCheck( stdoutRemove = True ):
        sys.exit("Execution Failed")
    logger.info("Genome index for BWA is generated.")

    
    
def starMapping(parsedArgs,reads,rlen,read_layout):
    logger = getLogger()
    overhang = rlen - 1
    
    star_dir = f"{parsedArgs.out_dir}/star_output"
    parsedArgs.mappingDir = star_dir
    if os.path.isdir(f"{star_dir}") == False:
        logger.info(f"Creating {star_dir} directory.")
        os.makedirs(f"{star_dir}")

    logger.info("Starting STAR mapping.")
    
    if parsedArgs.bacterial_mode is True:
        intron_len_max=1
    else:
        intron_len_max=0
        
    cmd_star_mapping = "STAR " \
f"--runThreadN {parsedArgs.threads} " \
f"--genomeDir {parsedArgs.starIndex} " \
f"--sjdbGTFfile {parsedArgs.out_dir}/concat.gtf " \
f"--sjdbOverhang {overhang} " \
f"--readFilesIn {reads} " \
"--readFilesCommand cat --outSAMtype BAM Unsorted " \
f"--outFileNamePrefix {star_dir}/concat_{rlen}_{read_layout}_ " \
f"--outFilterMismatchNmax {parsedArgs.outFilterMismatchNmax} " \
f"--outFilterMultimapNmax 10000 " \
f"--outFilterMismatchNoverReadLmax {parsedArgs.outFilterMismatchNoverReadLmax} " \
f"--alignIntronMax {intron_len_max} " \
f"--outTmpDir {parsedArgs.star_temp_dir}"

    # print(cmd_star_mapping)
    res = crossmapper.externalExec.execute(cmd_star_mapping,"STAR_mapping" , outDir = f"{star_dir}", overwrite = False)
    if not res.resCheck(clean = False):
        sys.exit("Execution Failed")
    
    logger.info("Mapping is finished. Started bam file sorting and indexing.")
    
    
    finalBamFile = f"{star_dir}/concat_{rlen}_{read_layout}_sorted.bam"
    tmpBamFile = f"{star_dir}/concat_{rlen}_{read_layout}_Aligned.out.bam"
    cmd_samtools_sort = "samtools sort " \
f"-@{parsedArgs.threads} " \
f"-o {finalBamFile} {tmpBamFile}"

    # print(cmd_samtools_sort)
    res  = crossmapper.externalExec.execute(cmd_samtools_sort,"samtools_sort" , outDir = f"{star_dir}")
    if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
        sys.exit("Execution Failed")
    
    logger.info("Sorting is finished. " +f"Final bam file writen to {finalBamFile}")
    
    
    try :
        logger.debug(f"Deleteing tmp sam file {tmpBamFile}")
        os.remove(tmpBamFile)
    except :
        logger.warning(f"Can not delete temporary bam file {tmpBamFile}")
    
    logger.info("Starting bam indexing.")
    cmd_samtools_index = f"samtools index {finalBamFile}"
    # print(cmd_samtools_index)
    res = crossmapper.externalExec.execute(cmd_samtools_index,"samtools_index" , outDir = f"{star_dir}")
    
    
    
    
    if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
        sys.exit("Execution Failed")
    logger.info("Bam Indexing is finished.")

        
        
    parsedArgs.mappingOutputFiles[rlen][read_layout] = finalBamFile
    
    
    
def bwaMapping(parsedArgs,reads,rlen,read_layout):
    logger = getLogger()
    logger.info("Starting mapping with BWA.")
#    cmd_bwa_mapping = "bwa mem " \
#    f"-t {parsedArgs.threads} concat {reads} -a | " \
#    f"samtools sort @{parsedArgs.threads} -o concat_{rlen}_{read_layout}_sorted.bam -"
#    print(cmd_bwa_mapping)
    
    bwa_dir = f"{parsedArgs.out_dir}/bwa_output"
    parsedArgs.mappingDir = bwa_dir
    if os.path.isdir(f"{bwa_dir}") == False:
        logger.info(f"Creating {bwa_dir} directory.")
        os.makedirs(f"{bwa_dir}")
    
    
    
    
    tmpSamFile = f"{bwa_dir}/concat_{rlen}_{read_layout}.sam"
    finalBamFile = f"{bwa_dir}/concat_{rlen}_{read_layout}_sorted.bam"
    cmd_bwa_mapping = f"bwa mem -a -t {parsedArgs.threads} -A {parsedArgs.match_score} -B {parsedArgs.mismatch_penalty} {parsedArgs.out_dir}/BWA_index/concat_BWA {reads}" 
    #f"samtools sort @{parsedArgs.threads} -o concat_{rlen}_{read_layout}_sorted.bam -
    res = crossmapper.externalExec.execute(cmd_bwa_mapping,"BWA_mapping" , 
                                 tmpSamFile,
                                 None,
                                 outDir = f"{bwa_dir}")
    if not res.resCheck():
        sys.exit("Execution Failed")
    ## TODO :: samtools view -bS
    res = crossmapper.externalExec.execute(f"samtools sort -@{parsedArgs.threads} -o {finalBamFile} {tmpSamFile}",
                                        "samtools" ,
                                        outDir = f"{bwa_dir}")
    if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
        sys.exit("Execution Failed")
    logger.info("Mapping is finished. " + f"Final bam file writen to {finalBamFile}")
    
    try :
        logger.debug(f"Deleteing tmp sam file {tmpSamFile}")
        os.remove(tmpSamFile)
    except :
        logger.warning(f"Can not delete temporary sam files {tmpSamFile}")
    
    
    logger.info("Starting Bam indexing.")
    cmd_samtools_index = f"samtools index {finalBamFile}"
    #print(cmd_samtools_index)
    res = crossmapper.externalExec.execute(cmd_samtools_index,
                                        "samtools_index",
                                        outDir = f"{bwa_dir}")
    if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
        sys.exit("Execution Failed")
    parsedArgs.mappingOutputFiles[rlen][read_layout] = finalBamFile
    logger.info("Bam Indexing is finished.")
    

    
def mapping(parsedArgs):
    
    if parsedArgs.simulation_type == "RNA":
        starIndex(parsedArgs)
    else:
        bwaIndex(parsedArgs)
    
    for rlen in parsedArgs.input_rlen:
        parsedArgs.mappingOutputFiles[rlen] = {}
        if parsedArgs.read_layout != "PE" : ## i.e SE or both
            # Single end case
            se_mapping = parsedArgs.simulationOutputFiles[rlen][0]
            if parsedArgs.simulation_type == "RNA":
                starMapping(parsedArgs,se_mapping,rlen,"SE")
            else:
                bwaMapping(parsedArgs,se_mapping,rlen,"SE")
        ####
        if parsedArgs.read_layout != "SE": ## i.e. PE or both
            pe_mapping = " ".join(parsedArgs.simulationOutputFiles[rlen])
            if parsedArgs.simulation_type == "RNA":
                starMapping(parsedArgs,pe_mapping,rlen,"PE")
            else:
                bwaMapping(parsedArgs,pe_mapping,rlen,"PE")
    
#    if parsedArgs.simulation_type == "RNA":
#        starIndex(parsedArgs)
#        for rlen in parsedArgs.input_rlen:
#            parsedArgs.mappingOutputFiles[rlen] = {}
#            se_mapping = os.path.join(f"{parsedArgs.out_dir}", f"concat_{rlen}_read1.fastq")
#            pe_mapping = os.path.join(f"{parsedArgs.out_dir}",f"concat_{rlen}_read1.fastq") + " " + os.path.join(f"{parsedArgs.out_dir}",f"concat_{rlen}_read2.fastq")
#            print(pe_mapping)
#            if parsedArgs.read_layout == "SE":
#                # TODO: cmd_remove_read2= "rm *_read2.fastq"
#                starMapping(parsedArgs,se_mapping,rlen,parsedArgs.read_layout)
#                
#            elif parsedArgs.read_layout == "PE":
#               starMapping(parsedArgs,pe_mapping,rlen,parsedArgs.read_layout)
#            else:
#                starMapping(parsedArgs,se_mapping,rlen,"SE")
#                starMapping(parsedArgs,pe_mapping,rlen,"PE")
#    else:
#        bwaIndex(parsedArgs)
#        for rlen in parsedArgs.input_rlen:
#            parsedArgs.mappingOutputFiles[rlen] = {}
#            #se_mapping = os.path.join(f"{parsedArgs.out_dir}", f"concat_{rlen}_read1.fastq")
#            se_mapping = parsedArgs.simulationOutputFiles[rlen][0]
#            #pe_mapping = os.path.join(f"{parsedArgs.out_dir}",f"concat_{rlen}_read1.fastq") + " " + os.path.join(f"{parsedArgs.out_dir}",f"concat_{rlen}_read2.fastq")
#            #pe_mapping = parsedArgs.simulationOutputFiles[rlen][0] + " " + parsedArgs.simulationOutputFiles[rlen][1]
#            if parsedArgs.read_layout == "SE":
#               #TODO: cmd_remove_read2= "rm *_read2.fastq"
#                bwaMapping(parsedArgs,se_mapping,rlen,parsedArgs.read_layout)
#            elif parsedArgs.read_layout == "PE":
#                bwaMapping(parsedArgs,
#                           " ".join(parsedArgs.simulationOutputFiles[rlen]),
#                           rlen,parsedArgs.read_layout)
#            else:
#                bwaMapping(parsedArgs,se_mapping,rlen,"SE")
#                bwaMapping(parsedArgs," ".join(parsedArgs.simulationOutputFiles[rlen]),rlen,"PE")
