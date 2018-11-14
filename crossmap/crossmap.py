# Imports
import sys
import os
import re
import sys
import argparse
import subprocess 
import math

from crossmap.helpers import getBaseName, setupLogger, getLogger , VerboseLevel
from crossmap.simulateReads import simulateData
from crossmap.mapping import prepareGenome
from crossmap.mapping import mapping



## temp allocation

_Main_Prog_Name = "crossmap.py" 
_Main_Prog_Desc = """
  --
  crossmap.py Software
"""

soft_version = "0.1"


standard_rlen = [25, 50, 75, 100, 125, 150, 300]

__DEBUG__ = True
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
    
    
    requirdSharedArgument =  shardParser.add_argument_group("Required Arguments")

    

    requirdSharedArgument.add_argument("-g", "--genomes",type=str, nargs=2, required=True,
    	help="Specify the genome files in fasta format. Enter genome names separated by whitespace. "
    	+ "\n NOTE: Keep the same order of listing for gtf/gff files")
    	
    requirdSharedArgument.add_argument("-a", "--annotations",type=str, nargs=2, required=True,
    	help="Specify the gtf/gff files. Enter the file names separated by whitespace. "
    	+ "NOTE: Keep the same order of listing as for genome files")
    	
    shardParser.add_argument("-t", "--threads", type=int, default = 1,
    	help = "Number of cores to be used for all multicore-supporting steps")
    	
    
    #### WGSIM parameters
    
    shardParser.add_argument("-e", "--error", type=float, default=0.02,
    	help = "Base error rate")
    	
    shardParser.add_argument("-d", "--outer_dist", type=int, default = 500,
    	help = 	"Outer distance between the two reads. For example, in case of 2x50 reads, d=300 and s=0 "
    			+ "the  mates will be 200 bp apart from each other.")
    
    shardParser.add_argument("-s", "--s_dev", type=int, default = 30,
    	help = "Standard deviation of outer distance.")
    	
    group = shardParser.add_mutually_exclusive_group(required=True)
    group.add_argument("-N", "--N_read", type = int, nargs=2,
        help = "The number of reads/read pairs to generate. This paremeter can not be used alongside with -C ")
    
    group.add_argument("-C", "--coverage", type = int, nargs=2,
        help = "Generate the number of reads that reaches the specified coverage. Coverage is calculated as:"
    		+ "C = N*rlen/L, " 
    		+ "where L is the length of the genome/transcriptome")
    	
    shardParser.add_argument("-rlay", "--read_layout", type = str, choices=["SE","PE","both"], default = "SE",
        help = "Specify the read configuration - single-end (SE), paired-end (PE), or both (both)." 
    		+ " If chosen 'both', the software will make separate analysis with each configuration")
    
    shardParser.add_argument("-rlen", "--read_length", type=str, default="50",
    	help = "Specify the read length. Choose from the possible read lengths available for Illumina machines:"
    		+ "25,50,75,100,125,150,300. The user can either enter a specific length, or specify a COMMA-SEPARATED (!)(no spaces are allowed between commas)"
    		+ "list of desired read lengths. In the latter case, the software will perform the analysis for all specified"
    		+ "values separatelly and will report mapping statistics in a form of a graph")
    		
    shardParser.add_argument("-r", "--mut_rate", type=float, default = 0.001,
    	help = "Mutation rate.")
    
    shardParser.add_argument("-R", "--indel_fraction", type=float, default = 0.015,
    	help = "Fraction of indels.")
    	
    shardParser.add_argument("-X", "--indel_extend", type=float, default = 0.3,
    	help = "Probability of an indel to be extended.")
    	
    shardParser.add_argument("-S", "--random_seed", type=int, default=(-1), 
    	help = "Seed for random generator.")
    
    shardParser.add_argument("-A", "--discard_ambig", type=float, default = 0.05,
    	help = "Disgard if the fraction of ambiguous bases is higher than this number.")
    	
    shardParser.add_argument("-hapl", "--haplotype_mode", action = "store_true", default = False,
    	help = "Haplotype mode. If specified, the haploid mutations will be simulated instead of diploid.")
    	
    
    shardParser.add_argument("-o", "--out_dir", default = "crossmap_out", type = str,
                       help = "Specify the output directory for crossmap output files.")
    
    
    #parsedArgs
    #parsedArgs = parser.parse_args()
    
    
    
    
    
    parser_DNA = subparsers.add_parser("DNA",help = "Simulate DNA data",formatter_class=argparse.ArgumentDefaultsHelpFormatter , parents=[shardParser] )
    
    dnaSharedGroup = parser_DNA.add_argument_group("Mapper Arguments","Arguments specific to BWA Mapper")
    dnaSharedGroup.add_argument("-specific_bwa", "--specific", type=int, default=5,
    	help = "Help message")
    	
    parser_RNA = subparsers.add_parser("RNA", help = "Simulate RNA data", formatter_class=argparse.ArgumentDefaultsHelpFormatter , parents=[shardParser])
    rnaSharedGroup = parser_RNA.add_argument_group("Mapper Arguments","Arguments specific to STAR Mapper")
    rnaSharedGroup.add_argument("-max_mismatch", "--outFilterMismatchNmax", type=int, default=10, metavar="Int",
    	help = "From STAR manual: "
    	+ " alignment will be output only if it has no more mismatches than this value")
    
    
    
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
        
        
    parsedArgs.fasta_names=[]
    if parsedArgs.simulation_type == "RNA":
        for i in range(0,len(parsedArgs.genomes)):
            transcriptome_name = getBaseName(parsedArgs.genomes[i]) + "_transcriptome%s"%(i+1) + ".fasta"
#            parsedArgs.fasta_names.append(os.path.abspath(transcriptome_name))
            parsedArgs.fasta_names.append( os.path.join(parsedArgs.out_dir,transcriptome_name))
            
    else:
        for i in range(0,len(parsedArgs.genomes)):
            parsedArgs.fasta_names.append(os.path.abspath(parsedArgs.genomes[i]))
    #print(parsedArgs.fasta_names)
    
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
    
    
    
    
#    if os.path.isdir("./%s"%(parsedArgs.out_dir)) == True:
#        print("%s already directory exists. Continuing."%(parsedArgs.out_dir))
#    else:
#        cmd_mkdir = "mkdir ./%s"%(parsedArgs.out_dir)
    

    
    
    ## other initilization 
    ## setting internal variable to parsedArgs object
    parsedArgs.isDebug = __DEBUG__
    parsedArgs.logPrefix = "crossmap.log"
    parsedArgs.logFile = os.path.join(parsedArgs.out_dir, parsedArgs.logPrefix)
    parsedArgs.verbose = VerboseLevel.All
    
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
    print("#" * 25)
    print(f"Simulation Type {parsedArgs.simulation_type}")
    print("Input Genomes Files :")
    indLevel = "\t"
    for inFile in parsedArgs.genomes:
        print(indLevel + "* " +  os.path.abspath(inFile))
    
    print("Species Prefix/Name and Id:")
    for speciesPrefix in parsedArgs.speciesPrefix:
        print(f"{indLevel} * {speciesPrefix}:{parsedArgs.speciesIds[speciesPrefix]}")
    if parsedArgs.simulation_type == "RNA" :
        print("Input Annotations Files :")
        indLevel = "\t"
        for inFile in parsedArgs.annotations:
            print(indLevel + "* " + os.path.abspath(inFile))
        print("Generated Transciptome Files :")
        indLevel = "\t"
        for inFile in parsedArgs.fasta_names:
            print(indLevel + "* " + os.path.abspath(inFile))  
    print(f"Read Length : {parsedArgs.input_rlen}")
    ## TODO :: complete the rest here
    

    return


## Main Script Entry Point
def crossmapMain():
    print("crossmapMain" , __package__)
    #print("Corssmap Test Run ...")
    #print(sys.argv)
    parser =  createArgumentParser()
    
    ## parse and check argument , also TODO :: may it would agood idea to prepare all parapmeters here if needed
    parsedArgs = parseArgument(parser)

    simulateData(parsedArgs)
    
    
    prepareGenome(parsedArgs)
    mapping(parsedArgs)

    return
