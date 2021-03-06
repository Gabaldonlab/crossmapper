#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 08:04:36 2019

@author: ahmed
"""
## add local module dir for easy import 

import sys
import os 

## get the module folder from the test script when executing it as script ...
crossmapMuduleDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__) ) +"/.." )
## or set it manually if running in spyder or ipython console 
## crossmapMuduleDir = "/home/ahmed/git/crossmap"
sys.path.append(crossmapMuduleDir)

import crossmapper 
import crossmapper.mapping


#%%
## to reload module with new changes
import importlib
importlib.reload(crossmapper)


#%% just test parser creatation
testCaseIndex = 3
testCases = [
        "crossmap.py RNA -h",
        ## 1
        "crossmap.py RNA -gb -rc -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta -a /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.93.gtf /data/bio/projects/simulation/input/test_11/Sbay_rebuilt.gtf -N 100000 100000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase11RNA -t 10",
        ## 2
        "crossmap.py DNA -gb -rc -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta  -N 1000 1000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase11DNA -t 10",
       
        ## 3 - test STAR template
        "crossmap.py RNA -gb -rc  --mapper-template /home/ahmed/git/crossmapper/config_examples/star.yaml -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta -a /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.93.gtf /data/bio/projects/simulation/input/test_11/Sbay_rebuilt.gtf -N 100000 100000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase_STAR_RNA -t 10",

        ## 4 - test BWA template
        "crossmap.py DNA -gb -rc --mapper-template /home/ahmed/git/crossmapper/config_examples/bwa.yaml  -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta  -N 100000 100000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase_BWA_DNA -t 10",
        ## 5 - test Bowtie2 template
        "crossmap.py DNA -gb -rc --mapper-template /home/ahmed/git/crossmapper/config_examples/bowtie2.yaml  -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta  -N 100000 100000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase_BOWTIE2_DNA -t 10",

        ## 6 - test HISAT2 template
        "crossmap.py RNA -gb -rc  --mapper-template /home/ahmed/git/crossmapper/config_examples/hisat2.yaml -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta -a /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.93.gtf /data/bio/projects/simulation/input/test_11/Sbay_rebuilt.gtf -N 100000 100000 -gn SC SU -rlen 50,100 -rlay both  -o /data/bio/projects/simulation/test/testcase_HISAT2_RNA -t 10",
        
        ## 7  - test STAR template
        
        "crossmap.py RNA -gb -rc --mapper-template /home/ghovhannisyan/users/tg/hhovhannisyan/crossmap_gabaldon_lab/crossmapper/config_examples/star.yaml -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -N 100000 100000 -rlen 50,100 -gn calb cpar -rlay both -o ./calb_cpar_star",
        
        
        ## 8 -  test BWA template
        
        "crossmap.py DNA -gb -rc --mapper-template /home/ghovhannisyan/users/tg/hhovhannisyan/crossmap_gabaldon_lab/crossmapper/config_examples/bwa.yaml -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 100000 100000 -rlen 50,100 -gn calb cpar -rlay both -o ./calb_cpar_bwa",
      
        ## 9 - test HISAT2 template
        
        "crossmap.py RNA -gb -rc --mapper-template /home/ghovhannisyan/users/tg/hhovhannisyan/crossmap_gabaldon_lab/crossmapper/config_examples/hisat2.yaml -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -N 100000 100000 -rlen 50,100 -gn calb cpar -rlay both -o ./calb_cpar_hisat2",
      
        ## 10 -test Bowtie
        
        "crossmap.py DNA -gb -rc --mapper-template /home/ghovhannisyan/users/tg/hhovhannisyan/crossmap_gabaldon_lab/crossmapper/config_examples/bowtie2.yaml -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 100000 100000 -rlen 50,100 -gn calb cpar -rlay both -o ./calb_cpar_bowtie2",

        ## 11 -test NGM
        
        "crossmap.py DNA -gb -rc --mapper-template /home/ghovhannisyan/users/tg/hhovhannisyan/crossmap_gabaldon_lab/crossmapper/config_examples/ngm.yaml -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 100000 100000 -rlen 50,100 -gn calb cpar -rlay both -o ./calb_cpar_ngm"

        ]


def getArgv(i):
    #print(testCases[i].split(' '))
    args = [] 
    for r in testCases[i].split(' '):
        if r.strip() != '':
            args.append(r)
            
    
    return args

def testCreateArgumentParser(caseIndex = 0):
    import crossmapper.crossmapper
    sys.argv =  getArgv(caseIndex)
    parser = crossmapper.crossmapper.createArgumentParser()
    parsedArgs = crossmapper.crossmapper.parseArgument(parser)
    #print(parsedArgs)
    return parsedArgs



def testRunCrossmap():
    sys.argv =  getArgv(0)
    #print(sys.argv)
    crossmapper.crossmapMain()

def testSimulation(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Simulation Test")
    crossmapper.mapping.concatGeneomes(parsedArgs)
    crossmapper.simulateReads.simulateData(parsedArgs)
    crossmapper.crossmapper.printOutputFileInfo(parsedArgs,'wgsim')
    logger.info("Simulation Test END")


def testMapping(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Mapping Test")
    crossmapper.mapping.mapping(parsedArgs)
    crossmapper.crossmapper.printOutputFileInfo(parsedArgs,'mapping')
    logger.info("Mapping Test END")    

def testMapper(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Mapping Test")
    parsedArgs.mapper.run()
    crossmapper.crossmapper.printOutputFileInfo(parsedArgs,'mapping')
    logger.info("Mapping Test END")  


def testCountingStep(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Counting Test")

    res = crossmapper.countUtil.getReadCounters(parsedArgs)
    logger.info("Counting Test END")    
    return res
#testCreateArgumentParser()
#testRunCrossmap()

def renameChromosomes(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting renameChromosomes Test")
    res = crossmapper.simulateReads.renameChromosomes(parsedArgs)
    logger.info("renameChromosomes Test END")    
    return res

#%%
#sys.argv =  getArgv(5)

args = testCreateArgumentParser(testCaseIndex)
renameChromosomes(args)
testSimulation(args)
testMapper(args)
#testMapping(args)
res = testCountingStep(args)
#%%
# importlib.reload(crossmap)
#args = testCreateArgumentParser()
#res = testCountingStep(args)
# args.groupBarChart = False
# crossmap.reporting.createHTMLReport(res,args)
## write result




#%%
# importlib.reload(crossmap)
# res = testCountingStep(args)


#exit(0)
#import crossmap.crossmap
#sys.argv =  getArgv(1)
#parser = crossmap.crossmap.createArgumentParser()
#parsedArgs = crossmap.crossmap.parseArgument(parser)
#%% 
## testing parser and general CLI 
#sys.argv =  getArgv()
## run main crossmap 
#crossmap.crossmap()

