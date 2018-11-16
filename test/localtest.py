
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
        "crossmap.py -g ../human.fasta ../mouse.fasta RNA -a ../human.gff ../mouse.gtf -C 20 20",
        "crossmap.py RNA -g ../human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o myres"
        "crossmap.py"
        ]


def getArgv(i):
    print(testCases[i].split(' '))
    return testCases[i].split(' ')

#%% just test parser creatation
def testCreateArgumentParser():
    import crossmap.crossmap
    sys.argv =  getArgv(0)
    parser = crossmap.crossmap.createArgumentParser()
    parsedArgs = crossmap.crossmap.parseArgument(parser)
    #print(parsedArgs)



def testRunCrossmap():
    sys.argv =  getArgv(0)
    #print(sys.argv)
    crossmap.crossmapMain()

#testCreateArgumentParser()
#testRunCrossmap()
#%%
#exit(0)
import crossmap.crossmap
sys.argv =  getArgv(0)
parser = crossmap.crossmap.createArgumentParser()
parsedArgs = crossmap.crossmap.parseArgument(parser)
#%% 
## testing parser and general CLI 
#sys.argv =  getArgv()
## run main crossmap 
#crossmap.crossmap()

#%%

import crossmap.externalExec



cmd= "cat /home/ahmed/git/crossmap/test/mydir/C_alb_A_chr1_50_read2.fastq /home/ahmed/git/crossmap/test/mydir/CPAR_chr1_50_read2.fastq"
# cat /home/ahmed/git/crossmap/test/mydir/C_alb_A_chr1_50_read2.fastq /home/ahmed/git/crossmap/test/mydir/CPAR_chr1_50_read2.fastq > /home/ahmed/git/crossmap/test/mydir/concat_50_read2.fastq
print(cmd)
#%%
crossmap.externalExec.execute(cmd,"cat","/home/ahmed/git/crossmap/test/mydir/concat_50_read2.fastq", outDir="/home/ahmed/git/crossmap/test/mydir/")
#cmd=["ls","."]

#crossmap.externalExec.execute("ls .")

#%%
import subprocess
cmd = "gffread asda zxczxc"
cmdList = cmd.split(" ")
#process = None

outfile = open("/home/ahmed/git/crossmap/test/mydir/concat_50_read2.fastq", "w")
errorfile = open("/home/ahmed/git/crossmap/test/mydir/cat_stderr.txt", "w")
process = subprocess.run(cmdList, shell=False, stdout = outfile, stderr = errorfile, check = False)
outfile.close()
errorfile.close()

