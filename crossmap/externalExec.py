## code for running external tools


import os
import sys
import subprocess
from crossmap.helpers import getLogger 


class ExecRes(object):
    def __init__(self,cmd,process, softName="extr_cmd" , stdOutFile = None  , stdErrFile = None ):
        self.cmd = cmd
        self.process = process
        self.returnCode = 0
        if process != None :
            self.returnCode = process.returncode
        self.softName = softName
        self.stdOutFile = stdOutFile
        self.stdErrFile = stdErrFile
    def clean(self, stdoutRemove=False,stdErrRemove = True):
        logger = getLogger()
        if stdoutRemove:
            logger.debug(f"Deleteing {self.stdOutFile}")
            os.remove(self.stdOutFile)
        if stdErrRemove:
            logger.debug(f"Deleteing {self.stdErrFile}")
            os.remove(self.stdErrFile)
        return
        ## return trye if not problems an
    def resCheck(self, clean = True , stdoutRemove=False,stdErrRemove = True ):
        logger = getLogger()
        if self.returnCode != 0 :
            logger.error(f"Error running {self.cmd}")
            logger.error(f"See error log in {self.stdErrFile}") 
            return False
        if clean :
            self.clean(stdoutRemove=stdoutRemove,stdErrRemove = stdErrRemove)
        return True
        
        


def execute(cmd, softName="extr_cmd" , stdOutFile = None  , stdErrFile = None , outDir = "" , overwrite = True):
    mode = "w"
    if not overwrite:
        mode = 'a'
    logger = getLogger()
    logger.debug(f"Start Running {softName} CMD : {cmd}")
    ## remove white spaces
    cmd = cmd.strip()
    cmd_list = cmd.split(" ") ## split as list 
    cmd_list = list(filter(None,cmd_list))
    if stdOutFile == None :
        stdOutFile = os.path.join(outDir , softName +"_stdout.txt")
        
    if stdErrFile == None:
        stdErrFile = os.path.join(outDir , softName +"_stderr.txt")
    
    
    with open(stdOutFile, mode) as outfile, open(stdErrFile, mode) as errorfile:
        try:
            if not overwrite:
                outfile.write("#"*25 + "\n")
                errorfile.write("#"*25 + "\n")
                outfile.write(f"{cmd}\n" + "#"*25 + "\n")
                errorfile.write(f"{cmd}\n" + "#"*25 + "\n")
                outfile.flush()
                errorfile.flush()
            process = subprocess.run(cmd_list, shell=False, stdout = outfile, stderr = errorfile, check = False)
            logger.debug("Running CMD Return Code : " + str(process.returncode))

            if process.returncode != 0:
                logger.error(f"Can not excecute {softName} CMD : \"{cmd}\".")
            return ExecRes(cmd,process,softName, stdOutFile , stdErrFile) 
        except FileNotFoundError as no_file:

            logger.error("Error in execute CMD : NO SUCH FILE OR DIRECTORY", exc_info=True)
        except PermissionError as perm_denied:
            logger.error("PERMISSION DENIED", exc_info=True)

        except Exception as ex:
            logger.error("Error in execute CMD", exc_info=True)            #raise ex
            
            
        return ExecRes(cmd,None,softName, stdOutFile , stdErrFile) 

