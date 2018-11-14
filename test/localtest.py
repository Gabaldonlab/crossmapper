
## add local module dir for easy import 

import sys
sys.path.append("/Users/hrant/crossmap")
import crossmap 

#%%
## to reload module with new changes
import importlib
importlib.reload(crossmap)


#%% 
testCases = [
        "crossmap.py RNA -h",
        "crossmap.py RNA -g ../human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir2/mydir",
        "crossmap.py"
        ]


def getArgv(i):
    print(testCases[i].split(' '))
    return testCases[i].split(' ')

#%% just test parser creatation
# =============================================================================
# def testCreateArgumentParser():
#     sys.argv =  getArgv(0)
#     parser = crossmap.createArgumentParser()
#     parsedArgs = parser.parse_args()
#     print(parsedArgs)
# =============================================================================


def testRunCrossmap():
    sys.argv =  getArgv(1)
    #print(sys.argv)
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

import crossmap.externalExec



cmd=["cd",".."]

print(cmd)
crossmap.externalExec.execute("ls .. && ls . ","star1")
cmd=["ls","."]

#crossmap.externalExec.execute("ls .")

