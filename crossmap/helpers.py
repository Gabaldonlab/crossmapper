
import sys
import os
import logging
from enum import Enum
###################### HELPER Methods #########################################
def getBaseName(filename):
    if len(filename.split("."))>1:
        
        basename = '.'.join(os.path.basename(filename).split(".")[0:-1])
    else:
        sys.exit("Error: please check the extensions of your input files."
                 +"\nGenome files should have .fa, .fasta or .fsa extensions."
                 +"\nGenome annotations should have .gtf or .gff extensions.")
    return basename    



class VerboseLevel(Enum):
    No = 1
    All = 2

def setupLogger(args):
    logger = logging.getLogger("crossmap")
    if args.isDebug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    # create a file handler
    f_handler = logging.FileHandler(args.logFile,"w")
    #f_handler.setLevel(logging.INFO)
    
    # create a logging format %(name)s -
    formatter = logging.Formatter('[%(asctime)s] %(levelname)-10s : %(message)s' ,  '%Y-%m-%d %H:%M:%S')
    f_handler.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.handlers = []
    logger.addHandler(f_handler)
    
    ##
    if args.verbose ==  VerboseLevel.All :
        c_handler = logging.StreamHandler()
        c_handler.setFormatter(formatter)
        logger.addHandler(c_handler)
    logger.propagate = False
    logger.debug("Debug is ON")
def getLogger():
    return logging.getLogger("crossmap")
