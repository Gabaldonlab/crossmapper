
## add local module dir for easy import 

import sys
import os 

## get the module folder from the test script when executing it as script ...
crossmapMuduleDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__) ) +"/.." )
## or set it manually if running in spyder or ipython console 
## crossmapMuduleDir = "/home/ahmed/git/crossmap"
sys.path.append(crossmapMuduleDir)
import crossmap 


#%%
## to reload module with new changes
import importlib
importlib.reload(crossmap)


#%% 
testCases = [
        "crossmap.py DNA  -h",
        "crossmap.py RNA -g ../human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o myres"
        ]
def getArgv(i = 0):
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
    crossmap.crossmapMain()

#testCreateArgumentParser()
#testRunCrossmap()
#%%
#exit(0)
import crossmap.crossmap
sys.argv =  getArgv(1)
parser = crossmap.crossmap.createArgumentParser()
parsedArgs = crossmap.crossmap.parseArgument(parser)
#%% 
## testing parser and general CLI 
#sys.argv =  getArgv()
## run main crossmap 
#crossmap.crossmap()

#%%