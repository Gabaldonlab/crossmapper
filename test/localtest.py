
## add local module dir for easy import 

import sys
sys.path.append("/home/ahmed/git/crossmap")
import crossmap 

#%%
## to reload module with new changes
import importlib
importlib.reload(crossmap)


#%% 
testCases = [
        "crossmap.py DNA  -h",
        "crossmap.py DNA -g ../human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir"
        ]
def getArgv(i = 0):
    return testCases[i].split(' ')

#%% just test parser creatation
def testCreateArgumentParser():
    sys.argv =  getArgv(0)
    parser = crossmap.createArgumentParser()
    parsedArgs = parser.parse_args()
    print(parsedArgs)


def testRunCrossmap():
    sys.argv =  getArgv(1)
    crossmap.crossmapMain()


#testCreateArgumentParser()
testRunCrossmap()

#exit(0)

#%% 
## testing parser and general CLI 
#sys.argv =  getArgv()
## run main crossmap 
#crossmap.crossmap()

#%%