
## add local module dir for easy import 

import sys
sys.path.append("/home/ahmed/git/crossmap/crossmap")
import crossmap

#%%
## to reload module with new changes
import importlib
importlib.reload(crossmap)


#%% 
testCases = [
        "crossmap.py DNA  -h",
        "crossmap.py RNA -g ../human.fasta ../../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir"
        ]
def getArgv(i = 0):
    return testCases[i].split(' ')

#%% just test parser creatation
def testCreateArgumentParser():
    sys.argv =  getArgv(1)
    parser = crossmap.createArgumentParser()
    parsedArgs = parser.parse_args()
    print(parsedArgs)


#%tb
#testCreateArgumentParser()
sys.argv =  getArgv(1)
parser = crossmap.createArgumentParser()
#parser = crossmap.ca()
try :    
    parsedArgs = crossmap.parseArgument(parser)
except Exception as ex:
    print("exit" , ex)
print(parsedArgs)

#exit(0)

#%% 
## testing parser and general CLI 
#sys.argv =  getArgv()
## run main crossmap 
#crossmap.crossmap()

#%%