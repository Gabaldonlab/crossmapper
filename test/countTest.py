
## add local module dir for easy import 

import sys
import os 

## get the module folder from the test script when executing it as script ...
crossmapMuduleDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__) ) +"/.." )
## or set it manually if running in spyder or ipython console 
## crossmapMuduleDir = "/home/ahmed/git/crossmap"
sys.path.append(crossmapMuduleDir)

import crossmap 
import crossmap.mapping


#%%
## to reload module with new changes
import importlib
importlib.reload(crossmap)


#%% 
testCases = [
        "crossmap.py DNA  -h",
        "crossmap.py RNA -g ./human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o myres",
        "crossmap.py DNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir -rlen 50,100", ## fixed coverage -- ~ same size,
        "crossmap.py RNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir -rlen 50,100", ## fail missing input files
        "crossmap.py RNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ./testFiles/C_alb_A_chr1.gff ./testFiles/CPAR_chr1.gff -C 20 20 -o mydir -rlen 50,100", ## 
        #"crossmap.py DNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -o mydir -rlen 50,100" ## fixed coverage -- ~ same size

        ]


def getArgv(i):
    #print(testCases[i].split(' '))
    return testCases[i].split(' ')

#%% just test parser creatation
def testCreateArgumentParser():
    import crossmap.crossmap
    sys.argv =  getArgv(1)
    parser = crossmap.crossmap.createArgumentParser()
    parsedArgs = crossmap.crossmap.parseArgument(parser)
    #print(parsedArgs)



def testRunCrossmap():
    sys.argv =  getArgv(1)
    #print(sys.argv)
    crossmap.crossmapMain()

def testSimulation():
    sys.argv =  getArgv(4)
    parser = crossmap.crossmap.createArgumentParser()
    parsedArgs = crossmap.crossmap.parseArgument(parser)
    logger = crossmap.helpers.getLogger()
    logger.info("Starting Simulation Test")
    crossmap.mapping.concatGeneomes(parsedArgs)
    crossmap.simulateReads.simulateData(parsedArgs)
    logger.info("Simulation Test END")

#testCreateArgumentParser()
#testRunCrossmap()
testSimulation()
#%%
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

