
import sys
import os
###################### HELPER Methods #########################################
def getBaseName(filename):
    if len(filename.split("."))>1:
        
        basename = '.'.join(os.path.basename(filename).split(".")[0:-1])
    else:
        sys.exit("Error: please check the extensions of your input files."
                 +"\nGenome files should have .fa, .fasta or .fsa extensions."
                 +"\nGenome annotations should have .gtf or .gff extensions.")
    return basename    
