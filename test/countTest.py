
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
        "crossmap.py RNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -N 1000 1000 -o testcase1 -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 1000 1000 -o testcase1 -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase2/CDUBL.fasta -N 1000 4000 -o testcase2DNA -rlen 25,50", ## 
        "crossmap.py RNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase2/CDUBL.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase2/CDUBL.gff -N 1000 4000 -o testcase2RNA -rlen 25,50", ## 
        "crossmap.py RNA -rlay both -g ./testFiles/testcase2/C_alb_A.fa ./testFiles/testcase2/CDUBL.fasta -a ./testFiles/testcase2/C_alb_A.gff ./testFiles/testcase2/CDUBL.gff -C 2 2  -o /data/bio/projects/simulation/test/testcase3RNA -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase3/C_alb_A_chr1.fasta ./testFiles/testcase3/CPAR_chr1.fasta ./testFiles/testcase3/dubl.fasta  -C 2 2 2 -o /data/bio/projects/simulation/test/testcase4DNA -rlen 50,100", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase3/C_alb_A_chr1.fasta ./testFiles/testcase3/CPAR_chr1.fasta ./testFiles/testcase3/dubl.fasta  -N 100000 100000 100000 -o /data/bio/projects/simulation/test/testcase5DNA -rlen 50,100", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase4/C_alb_A_chr1.fasta ./testFiles/testcase4/CPAR_chr1.fasta ./testFiles/testcase4/dubl.fasta  -N 100000 100000 100000 -o /data/bio/projects/simulation/test/testcase6DNA -rlen 50,100,125", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -N 1000000 1000000 1000000 -o /data/bio/projects/simulation/test/testcase7DNA -rlen 50,100,125", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -N 3000000 1000000 500000  -o /data/bio/projects/simulation/test/testcase8DNA -rlen 50,125", ## 
        "crossmap.py DNA  -gb -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -C 3 1 1  -o /data/bio/projects/simulation/test/testcase9DNA -rlen 125", ## 

        #"crossmap.py DNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -o mydir -rlen 50,100" ## fixed coverage -- ~ same size

        ]


def getArgv(i):
    #print(testCases[i].split(' '))
    args = [] 
    for r in testCases[i].split(' '):
        if r.strip() != '':
            args.append(r)
            
    
    return args

#%% just test parser creatation
def testCreateArgumentParser():
    import crossmap.crossmap
    sys.argv =  getArgv(9)
    parser = crossmap.crossmap.createArgumentParser()
    parsedArgs = crossmap.crossmap.parseArgument(parser)
    crossmap.mapping.concatGeneomes(parsedArgs)
    #print(parsedArgs)
    return parsedArgs



def testRunCrossmap():
    sys.argv =  getArgv(1)
    #print(sys.argv)
    crossmap.crossmapMain()

def testSimulation(parsedArgs):
    logger = crossmap.helpers.getLogger()
    logger.info("Starting Simulation Test")
    crossmap.simulateReads.simulateData(parsedArgs)
    crossmap.crossmap.printOutputFileInfo(parsedArgs,'wgsim')
    logger.info("Simulation Test END")


def testMapping(parsedArgs):
    logger = crossmap.helpers.getLogger()
    logger.info("Starting Mapping Test")
    crossmap.mapping.mapping(parsedArgs)
    crossmap.crossmap.printOutputFileInfo(parsedArgs,'mapping')
    logger.info("Mapping Test END")    

def testCountingStep(parsedArgs):
    logger = crossmap.helpers.getLogger()
    logger.info("Starting Counting Test")

    res = crossmap.countUtil.getReadCounters(parsedArgs)
    logger.info("Counting Test END")    
    return res
#testCreateArgumentParser()
#testRunCrossmap()
    
    
#sys.argv =  getArgv(5)

args = testCreateArgumentParser()
testSimulation(args)
testMapping(args)
res = testCountingStep(args)
#%%
#importlib.reload(crossmap)
#args = testCreateArgumentParser()
#res = testCountingStep(args)
# args.groupBarChart = False
#crossmap.reporting.createHTMLReport(res,args)
## write result

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

