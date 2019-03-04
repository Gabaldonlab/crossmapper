# Imports
import sys
import os
import re
import sys
import argparse
import subprocess 
import math

from crossmap.helpers import getBaseName, setupLogger, getLogger , VerboseLevel
from crossmap.simulateReads import simulateData, renameChromosomes
from crossmap.mapping import concatGeneomes
from crossmap.mapping import mapping
from crossmap.countUtil import getReadCounters


## temp allocation

_Main_Prog_Name = "crossmapper" 
_Main_Prog_Desc = """
  --
  Crossmapper Software
"""

soft_version = "0.1"


standard_rlen = [50, 75, 100, 125, 150, 300]

__DEBUG__ = False
###############################################################################



###################### Argument Parsing and checking ##########################
def createArgumentParser():
#TODO: Now the subparses only can appear as the last argumnet, need to fix it

###Create top level parser

#formatter_class=argparse.RawTextHelpFormatter - with this the text in help is free of formating
#formatter_class=argparse.ArgumentDefaultsHelpFormatter - this adds default values in help. Should be in all subparsers.

    mainParser = argparse.ArgumentParser(prog = _Main_Prog_Name,description = _Main_Prog_Desc , formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    
    mainParser.add_argument("-v", "--version", action = "version", \
        version = "%(prog)s \"v" + soft_version + "\"")
    
    ## subparser for DNA or RNA running mode
    subparsers = mainParser.add_subparsers(help="Simulation type. Choose to simulate either DNA or RNA data",title = "SimulationType" ,  dest = "simulation_type")
    subparsers.required = True
    
    
    
    shardParser = argparse.ArgumentParser(add_help=False)
    
    
    requirdSharedArgument = shardParser.add_argument_group("Required Arguments")

    requirdSharedArgument.add_argument("-g", "--genomes",type=str, nargs="+", required=True, metavar = "fasta",
    	help="Specify the genome files in fasta format. Enter genome names separated by whitespace. "
    	+ "\n NOTE: Keep the same order of listing for gtf/gff files")
    
    shardParser.add_argument("-t", "--threads", type=int, default = 1, metavar = "int",
    	help = "Number of cores to be used for all multicore-supporting steps")
    	
    
    #### WGSIM parameters
    
    shardParser.add_argument("-e", "--error", type=float, default=0.02, metavar = "float",
    	help = "Base error rate")
    	
    shardParser.add_argument("-d", "--outer_dist", type=int, default = 500, metavar = "int",
    	help = 	"Outer distance between the two reads. For example, in case of 2x50 reads, d=300 and s=0 "
    			+ "the  mates will be 200 bp apart from each other.")
    
    shardParser.add_argument("-s", "--s_dev", type=int, default = 30, metavar = "int",
    	help = "Standard deviation of outer distance.")
    	
    group = shardParser.add_mutually_exclusive_group(required=True)
    group.add_argument("-N", "--N_read", type = int, nargs="+",metavar = "int", 
        help = "The number of reads/read pairs to generate. This parameter can not be used alongside with -C ")
    
    group.add_argument("-C", "--coverage", type = float, nargs="+",metavar="float/int",
        help = "Generate the number of reads that reaches the specified coverage. Coverage is calculated as:"
    		+ "C = N*rlen/L, " 
    		+ "where L is the length of the genome/transcriptome")
    	
    shardParser.add_argument("-rlay", "--read_layout", type = str, choices=["SE","PE","both"], default = "SE",
        help = "Specify the read configuration - single-end (SE), paired-end (PE), or both (both)." 
    		+ " If chosen 'both', the software will make separate analysis with each configuration")
    
    shardParser.add_argument("-rlen", "--read_length", type=str, default="50", metavar = "int",
    	help = "Specify the read length. Choose from the possible read lengths available for Illumina machines: "
    		+ "50,75,100,125,150,300. The user can either enter a specific length, or specify a (!) COMMA-SEPARATED (no spaces are allowed between commas) "
    		+ "list of desired read lengths. In the latter case, the software will perform the analysis for all specified "
    		+ "values separatelly and will report mapping statistics in a form of a graph")
    		
    shardParser.add_argument("-r", "--mut_rate", type=float, default = 0.001,metavar="float",
    	help = "Mutation rate.")
    
    shardParser.add_argument("-R", "--indel_fraction", type=float, default = 0.015,metavar="float",
    	help = "Fraction of indels.")
    	
    shardParser.add_argument("-X", "--indel_extend", type=float, default = 0.3,metavar="float",
    	help = "Probability of an indel to be extended.")
    	
    shardParser.add_argument("-S", "--random_seed", type=int, default=(-1), metavar="int",
    	help = "Seed for random generator.")
    
    shardParser.add_argument("-AMB", "--discard_ambig", type=float, default = 0.05,metavar="float",
    	help = "Disgard if the fraction of ambiguous bases is higher than this number.")
    	
    shardParser.add_argument("-hapl", "--haplotype_mode", action = "store_true", default = False,
    	help = "Haplotype mode. If specified, the haploid mutations will be simulated instead of diploid.")
    	
    
    shardParser.add_argument("-o", "--out_dir", default = "crossmap_out", type = str, metavar = "PATH",
                       help = "Specify the output directory for crossmap output files.")


    ## TODO ::
    ## add args.graphPercent option
    ## groupBarChart
    ## add args.groupBarChart option
    shardParser.add_argument("-gb", "--groupBarChart",  action = "store_true", default = False,
                       help = "Use a grouped bar chart in the output report instead of individual bar chart.")
    
    shardParser.add_argument("-rc", "--reportCrossmapped",  action = "store_true", default = False,
                       help = "Report all cross mapped reads into csv file.")
    
    shardParser.add_argument("-gn", "--genome_names", type=str, nargs="+", metavar = "name",
                                      help="Specify names of the genomes. The names will appear in the report file." )

    parser_DNA = subparsers.add_parser("DNA",help = "Simulate DNA data",formatter_class=argparse.ArgumentDefaultsHelpFormatter , parents=[shardParser] )
    
    dnaSharedGroup = parser_DNA.add_argument_group("Mapper Arguments","Arguments specific to BWA Mapper")
    dnaSharedGroup.add_argument("-k", "--seed_length", type=int, default=19,metavar="int",
    	help = "From BWA manual: "
        +" Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20")
    
    dnaSharedGroup.add_argument("-A", "--match_score", type=int, default=1, metavar="int",
                                help = "From BWA manual: "
                                + " Matching score.")
    dnaSharedGroup.add_argument("-B", "--mismatch_penalty", type=int, default=4, metavar="int",
                                help = "From BWA manual: "
                                + " Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}")
    	
    parser_RNA = subparsers.add_parser("RNA", help = "Simulate RNA data", formatter_class=argparse.ArgumentDefaultsHelpFormatter , parents=[shardParser])
    rnaSharedGroup = parser_RNA.add_argument_group("Mapper and annotation Arguments","Arguments specific to STAR Mapper")
    
    rnaSharedGroup.add_argument("-max_mismatch_per_len", "--outFilterMismatchNoverReadLmax", type=float, default=0.04, metavar = "float",
                                help = "From STAR manual: "
                                  +"alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value: for 2x100b, max number of mismatches is 0.04*200=8 for the paired read.")
    rnaSharedGroup.add_argument("-bact_mode", "--bacterial_mode", action = "store_true", default=False,
                                help = "This option prohibits spliced alignments for STAR and it can be used for mapping bacterial data." )

    
    
    rnaSharedGroup.add_argument("-max_mismatch", "--outFilterMismatchNmax", type=int, default=10, metavar="int",
    	help = "From STAR manual: "
    	+ " alignment will be output only if it has no more mismatches than this value")
    
    rnaSharedGroup.add_argument("-a", "--annotations",type=str, nargs="+", required=True,metavar = "gtf",
    	help="Specify the gtf/gff files. Enter the file names separated by whitespace. "
    	+ "NOTE: Keep the same order of listing as for genome files")
    

    rnaSharedGroup.add_argument("-star_tmp", "--star_temp_dir", default = "./TMPs",  type=str, metavar="PATH",
    	help = "Specify a full path to a local temprorary directory, where all intermediate files of STAR will be written. "
    	+ " This option can be used when running Crossmaper from a local machine in a server or cluster with SAMBA connection.")

    
    
    return mainParser


def parseArgument(argumentParser):
    parsedArgs = argumentParser.parse_args()

    ## setup absole path for dir
    parsedArgs.out_dir = os.path.abspath(parsedArgs.out_dir) 
    if os.path.isdir(parsedArgs.out_dir) != True:
        ## TODO :: create the folder here
        # cmd_mkdir = "mkdir ./%s"%(parsedArgs.out_dir)
        ## try and handle execption here
        os.makedirs(parsedArgs.out_dir)


    for i in range(0,len(parsedArgs.genomes)):
        if os.path.exists(parsedArgs.genomes[i]):
            if not os.path.getsize(parsedArgs.genomes[i]) > 0:
                sys.exit(f"Error: {parsedArgs.genomes[i]} file is empty! Please provide a valid file.")
        else:
            sys.exit(f"Error: {parsedArgs.genomes[i]} file does not exist! Please provide a valid file.")
                
############### checking input
    if len(parsedArgs.genomes) <= 1 :
            sys.exit(f"Error: Number of provided input genomes must be at least 2.")
        
    if parsedArgs.simulation_type == "RNA":
        if len(parsedArgs.genomes) != len(parsedArgs.annotations):
            sys.exit(f"Error: Number of provided input genomes files does not match number of input annotations files.")

    if parsedArgs.coverage is not None:
        if len(parsedArgs.coverage) == 1 :
            for ic in range(1,len(parsedArgs.genomes)):
                parsedArgs.coverage.append(parsedArgs.coverage[0])
        if len(parsedArgs.genomes) > len(parsedArgs.coverage):
            sys.exit(f"Error: Provided Coverage (--coverage) options do not match the input genomes files. You should provide coverage for each input fasta file or just one coverage for all of them.")
    elif parsedArgs.N_read is not None:
        if len(parsedArgs.N_read) == 1 :
            for ic in range(1,len(parsedArgs.genomes)):
                parsedArgs.N_read.append(parsedArgs.N_read[0])
        elif len(parsedArgs.genomes) > len(parsedArgs.N_read):
            sys.exit(f"Error: Provided  number of reads/read pairs to generate (--N_read) options do not match the input genomes files. You should provide one for each input fasta file or just one for all of them.")
       

    ### for renaming chr names
    parsedArgs.chr_rename_fasta=[]
    for i in range(0,len(parsedArgs.genomes)):
        fasta_chr_rename=getBaseName(parsedArgs.genomes[i])+"_rename"+".fasta"
        parsedArgs.chr_rename_fasta.append(os.path.abspath(parsedArgs.out_dir)+"/"+fasta_chr_rename)
    #print(parsedArgs.chr_rename_fasta)
    
    
    ### for renaming chr names in gff
    if parsedArgs.simulation_type == "RNA":
        parsedArgs.chr_rename_gff=[]
        for i in range(0,len(parsedArgs.annotations)):
            if parsedArgs.annotations[i][-3:] == "gtf":
                gff_chr_rename=getBaseName(parsedArgs.annotations[i])+"_rename"+".gtf"
                parsedArgs.chr_rename_gff.append(os.path.abspath(parsedArgs.out_dir)+"/"+gff_chr_rename)
            elif parsedArgs.annotations[i][-3:] == "gff":
                gff_chr_rename=getBaseName(parsedArgs.annotations[i])+"_rename"+".gff"
                parsedArgs.chr_rename_gff.append(os.path.abspath(parsedArgs.out_dir)+"/"+gff_chr_rename)
            
        #print(parsedArgs.chr_rename_gff)



    parsedArgs.fasta_names=[]
    if parsedArgs.simulation_type == "RNA":
        for i in range(0,len(parsedArgs.chr_rename_fasta)):
            transcriptome_name = getBaseName(parsedArgs.chr_rename_fasta[i]) + "_transcriptome%s"%(i+1) + ".fasta"
#            parsedArgs.fasta_names.append(os.path.abspath(transcriptome_name))
            parsedArgs.fasta_names.append(os.path.join(parsedArgs.out_dir,transcriptome_name))

        if len(parsedArgs.annotations)>0:
            for i in range(0,len(parsedArgs.annotations)):
                if os.path.exists(parsedArgs.annotations[i]):
                    if not os.path.getsize(parsedArgs.annotations[i]) > 0:
                        sys.exit(f"Error: {parsedArgs.annotations[i]} file is empty! Please provide a valid file.")
                else:
                    sys.exit(f"Error: {parsedArgs.annotations[i]} file does not exist! Please provide a valid file.")
            
    else:
        for i in range(0,len(parsedArgs.chr_rename_fasta)):
            parsedArgs.fasta_names.append(os.path.abspath(parsedArgs.chr_rename_fasta[i]))
    #print(parsedArgs.fasta_names)
    
    
    ## check if not all values can be converted to int
    try:
        list(map(int,parsedArgs.read_length.split(",")))
    except Exception:
        sys.exit("There are strings or floats in read length values. Please use only standard read lengths!")
    
    ## convert list of strings to list of integers
    input_rlen=list(map(int,parsedArgs.read_length.split(",")))
    #print(input_rlen)
    
    ## check if there are duplicated lengths
    if not len(set(input_rlen)) == len(input_rlen):
        sys.exit("Error: read lengths shoud not be duplicated!")    
    
    ## check if any length is not standard
    for length in input_rlen:
        #print(length)
        if not length in standard_rlen:
            sys.exit("Error: input read length %s is not a standard Illumina read length."%(length)
                     + "\nPlease refer to our help page (crossmap -h) to find standard read lengths.")    
    
    
    parsedArgs.input_rlen = input_rlen
    
    


    
    
    ## other initilization 
    
    parsedArgs.simulationOutputFiles = {}
    parsedArgs.mappingOutputFiles = {}
    ## default concat geneome fasta file name
    parsedArgs.simDir = f"{parsedArgs.out_dir}"
    parsedArgs.mappingDir = f"{parsedArgs.out_dir}"
    parsedArgs.genomeConcatFasta = f"{parsedArgs.out_dir}/concat.fasta"
    parsedArgs.annotationsGTFConcat = f"{parsedArgs.out_dir}/concat.gtf"
    ## setting internal variable to parsedArgs object
    parsedArgs.isDebug = __DEBUG__
    parsedArgs.logPrefix = "crossmap.log"
    parsedArgs.logFile = os.path.join(parsedArgs.out_dir, parsedArgs.logPrefix)
    parsedArgs.verbose = VerboseLevel.All
    ## option to report crossmapp reads info files
    # parsedArgs.reportCrossmapped = True
    
    
    if parsedArgs.genome_names is not None:
        parsedArgs.speciesPrefix = parsedArgs.genome_names
    else:
        parsedArgs.speciesPrefix = None
    
    
    if parsedArgs.speciesPrefix == None :
        ## get basename from the genome file
        parsedArgs.speciesPrefix = []
        for i in range(0,len(parsedArgs.genomes)):
            genomePrefix = getBaseName(parsedArgs.genomes[i])
            parsedArgs.speciesPrefix.append(genomePrefix)
    
    ## create specied Ids dict
    parsedArgs.speciesIds = {}
    for i in range(0,len(parsedArgs.speciesPrefix)) :
         parsedArgs.speciesIds[parsedArgs.speciesPrefix[i]] = i
    
    
    
    if __DEBUG__:
        printArgs(parsedArgs)
    cmdLine = " ".join(sys.argv )
    setupLogger(parsedArgs)

    getLogger().info("Starting the program with  \"" + cmdLine + "\"")

    return parsedArgs





def printArgs(parsedArgs):
    allArgsStr = "All Args \n"
    allArgsStr += "#" * 25 + "\n"
    allArgsStr += f"Simulation Type {parsedArgs.simulation_type}" + "\n"
    allArgsStr += "Input Genomes Files :" + "\n"
    indLevel = "\t"
    for inFile in parsedArgs.genomes:
        allArgsStr += indLevel + "* " +  os.path.abspath(inFile) + "\n"
    
    allArgsStr += "Species Prefix/Name and Id:" + "\n"
    for speciesPrefix in parsedArgs.speciesPrefix:
        allArgsStr += f"{indLevel} * {speciesPrefix}:{parsedArgs.speciesIds[speciesPrefix]}" + "\n"
    if parsedArgs.simulation_type == "RNA" :
        allArgsStr +=  "Input Annotations Files :"+ "\n"
        indLevel = "\t"
        for inFile in parsedArgs.annotations:
            allArgsStr += indLevel + "* " + os.path.abspath(inFile)+ "\n"
        allArgsStr +="Generated Transciptome Files :"+ "\n"
        indLevel = "\t"
        for inFile in parsedArgs.fasta_names:
            allArgsStr += indLevel + "* " + os.path.abspath(inFile)+ "\n"
    allArgsStr += f"Read Length : {parsedArgs.input_rlen}"+ "\n"
    allArgsStr += f"Read Layout : {parsedArgs.read_layout}"+ "\n"
    allArgsStr += "#" * 25 + "\n"

    getLogger().debug(allArgsStr)
    ## TODO :: complete the rest here

    return


def printOutputFileInfo(parsedArgs,step = 'All'):
    logger = getLogger()
    if step == "wgsim" or step == 'All':
        allArgsStr ="\n"
        for rlen,files in parsedArgs.simulationOutputFiles.items():
            allArgsStr +=  f"\t\t * {rlen} : {files}\n"
        logger.debug("wgsim output files : "+ allArgsStr)
    if step == "mapping" or step == 'All':
        allArgsStr ="\n"
        for rlen,layout_files in parsedArgs.mappingOutputFiles.items():
            for layout,files in layout_files.items():
                allArgsStr +=  f"\t\t * {rlen} ({layout}) : {files}\n"
        logger.debug("mapping output files : "+ allArgsStr)
        
        


## Main Script Entry Point
def crossmapMain():
    print("crossmapMain" , __package__)
    #print("Corssmap Test Run ...")
    #print(sys.argv)
    parser =  createArgumentParser()
    
    ## parse and check argument , also TODO :: may it would agood idea to prepare all parapmeters here if needed
    parsedArgs = parseArgument(parser)
    

    
    ## 1 rename Chromosomes in fasta file and gtf file
    renameChromosomes(parsedArgs)
    
    ## 2 concat genomes files
    concatGeneomes(parsedArgs)
    
    ## 3 
    simulateData(parsedArgs)
    
    ## 4 
    mapping(parsedArgs)

    ## 5 saving files and report are handled inside 
    finalResult = getReadCounters(parsedArgs)

    return
