
## add local module dir for easy import 

import sys
import os 

## get the module folder from the test script when executing it as script ...
crossmapMuduleDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__) ) +"/.." )
## or set it manually if running in spyder or ipython console 
## crossmapMuduleDir = "/home/ahmed/git/crossmap"
sys.path.append(crossmapMuduleDir)

import crossmapper 
import crossmapper.mapping


#%%
## to reload module with new changes
import importlib
importlib.reload(crossmapper)


#%% 
testCases = [
        "crossmapper --version "
        "crossmapper DNA -g /data/bio/projects/simulation/input/bug/C_glabrata_CBS138_current_chromosomes.fasta /data/bio/projects/simulation/input/bug/C_parapsilosis_CDC317_current_chromosomes.fasta -t 8 -C 5 5 -rlay both -rlen 50 -gb -rc -gn CGLAB CPAR -o /data/bio/projects/simulation/test",
        "crossmap.py RNA -gb -rc -g  /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa /data/bio/projects/simulation/input/test_11/Sbay.ultrascaf_AND_unplaced.fasta -a /data/bio/projects/simulation/input/test_11/Saccharomyces_cerevisiae.R64-1-1.93.gtf /data/bio/projects/simulation/input/test_11/Sbay_rebuilt.gtf -N 1000 1000 -gn SC SU -rlen 50 -rlay PE  -o /data/bio/projects/simulation/test/testcase11RNA -t 10",
        "crossmap.py RNA -h",
        "crossmap.py RNA -g ./human.fasta ../mouse.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o myres",
        "crossmap.py DNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir -rlen 50,100", ## fixed coverage -- ~ same size,
        "crossmap.py RNA -g ./testFiles/C_alb_A_chr1.fasta ./testFiles/CPAR_chr1.fasta -a ../human.gff ../mouse.gtf -C 20 20 -o mydir -rlen 50,100", ## fail missing input files
        "crossmap.py RNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -N 1000 1000 -o testcase1 -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 1000 1000 -o testcase1 -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase2/CDUBL.fasta -N 1000 4000 -o testcase2DNA -rlen 25,50", ## 
        "crossmap.py RNA -rlay both -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase2/CDUBL.fasta -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase2/CDUBL.gff -N 1000000 1000000 -o /data/bio/projects/simulation/test/testcase2RNA -rlen 50,150", ## 
        "crossmap.py RNA -rlay both -g ./testFiles/testcase2/C_alb_A.fa ./testFiles/testcase2/CDUBL.fasta -a ./testFiles/testcase2/C_alb_A.gff ./testFiles/testcase2/CDUBL.gff -C 2 2  -o /data/bio/projects/simulation/test/testcase3RNA -rlen 25,50", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase3/C_alb_A_chr1.fasta ./testFiles/testcase3/CPAR_chr1.fasta ./testFiles/testcase3/dubl.fasta  -C 2 2 2 -o /data/bio/projects/simulation/test/testcase4DNA -rlen 50,100", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase3/C_alb_A_chr1.fasta ./testFiles/testcase3/CPAR_chr1.fasta ./testFiles/testcase3/dubl.fasta  -N 100000 100000 100000 -o /data/bio/projects/simulation/test/testcase5DNA -rlen 50,100", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase4/C_alb_A_chr1.fasta ./testFiles/testcase4/CPAR_chr1.fasta ./testFiles/testcase4/dubl.fasta  -N 100000 100000 100000 -o /data/bio/projects/simulation/test/testcase6DNA -rlen 50,125", ## 
        "crossmap.py DNA -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -N 1000000 1000000 1000000 -o /data/bio/projects/simulation/test/testcase7DNA -rlen 50,100,125", ## 
        "crossmap.py DNA  -rc -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -N 3000000 2000000 200000  -o cross_report -rlen 50,125", ## 
        "crossmap.py DNA  -gb -rlay both -g ./testFiles/testcase5/C_alb_A_chr1.fasta ./testFiles/testcase5/CPAR_chr1.fasta ./testFiles/testcase5/dubl.fasta  -C 3 1 1  -o /data/bio/projects/simulation/test/testcase9DNA -rlen 125", ## 
        "crossmap.py RNA -r 0.01 -max_mismatch 5 -gb -rlay both -g ./testFiles/testcase3/C_alb_A.fasta ./testFiles/testcase3/C_alb_B.fasta -a ./testFiles/testcase3/C_alb_A.gff ./testFiles/testcase3/C_alb_B.gff -N 1000000 1000000 -o /home/ghovhannisyan/users/tg/hhovhannisyan/crossmaping_tool/master_script/test/test_CALB_A_vs_B_mism_10 -rlen 50,75,100,125,150 -star_tmp /home/ghovhannisyan/TMP/TMP3 -t 10",
        "crossmap.py DNA -r 0.01 -gb -rlay both -g ./testFiles/test_homo_mus_ceano_dros/homo_rename.fasta ./testFiles/test_homo_mus_ceano_dros/mus_rename.fasta ./testFiles/test_homo_mus_ceano_dros/dros_rename.fasta ./testFiles/test_homo_mus_ceano_dros/caeno_rename.fasta -N 1000000 1000000 1000000 1000000 -o /home/ghovhannisyan/users/tg/hhovhannisyan/crossmaping_tool/master_script/test/4_species_fixed -rlen 50,75,100,125,150 -t 8",
        "crossmap.py DNA -r 0.01 -gb -rlay SE -g /data/bio/data/sim/testcase10/homo_.fasta /data/bio/data/sim/testcase10/mus_.fasta /data/bio/data/sim/testcase10/dros_.fasta /data/bio/data/sim/testcase10/caeno_.fasta -N 1000000 500000 500000 500000 -o /data/bio/projects/simulation/test/4_species -rlen 50 -t 8",
        "crossmap.py RNA -r 0.01 -max_mismatch 5 -gb -rlay SE -g ./testFiles/testcase3/C_alb_A.fasta ./testFiles/testcase3/C_alb_B.fasta -a ./testFiles/testcase3/C_alb_A.gff ./testFiles/testcase3/C_alb_B.gff -N 1000000 1000000 -o /data/bio/projects/simulation/test/testcase3RNA -star_tmp /data/bio/projects/simulation/test/testcase3RNA/TMP_  -rlen 50,75 -t 10",
        "crossmap.py RNA -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -gb -N 100000 100000 -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -o ./renaming_test -rlen 50,75 -star_tmp /home/ghovhannisyan/TMP22 -t 10", ## check renaming
        "crossmap.py DNA -gn 1 2 -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -N 100000 100000 -o ./renaming_test_DNA -rlen 50,75 -t 10", ## check renaming
        "crossmap.py RNA -g ./testFiles/testcase1/C_alb_A_chr1.fasta ./testFiles/testcase1/CPAR_chr1.fasta -gb -N 10000 10000 -a ./testFiles/testcase1/C_alb_A_chr1.gff ./testFiles/testcase1/CPAR_chr1.gff -o ./renaming_test -rlen 50,75 -star_tmp /home/ghovhannisyan/TMP22 -t 10", ## check renaming
        "crossmap.py DNA -r 0.01 -gb -rlay both -gn human mouse fly nematode -g ./testFiles/test_homo_mus_ceano_dros/homo.fasta ./testFiles/test_homo_mus_ceano_dros/mus.fasta ./testFiles/test_homo_mus_ceano_dros/dros.fasta ./testFiles/test_homo_mus_ceano_dros/caeno.fasta -N 2500000 2500000 2500000 2500000 -o /home/ghovhannisyan/users/tg/hhovhannisyan/crossmaping_tool/master_script/test/4_species_fixed_new_names_15Dec -rlen 50,75,100,125,150 -t 11",
        "crossmap.py RNA -gb -rlay both -rlen 50,75,100,125 -gn yeast human -g ./human_cglab/C_glabrata_CBS138_current_chromosomes.fasta ./human_cglab/Homo_sapiens.GRCh38.dna.primary_assembly.fa -a ./human_cglab/C_glabrata_CBS138_current_features.gff ./human_cglab/Homo_sapiens.GRCh38.94.gtf -N 5000000 5000000 -star_tmp /home/ghovhannisyan/TMP69 -t 7 -o /home/ghovhannisyan/users/tg/hhovhannisyan/crossmaping_tool/master_script/test/human_cglab_test",
        "crossmap.py RNA -gb -rlay both -rlen 50 -gn yeast human -g ./human_cglab/C_glabrata_CBS138_current_chromosomes.fasta ./human_cglab/human_chr21.fasta -a ./human_cglab/C_glabrata_CBS138_current_features.gff ./human_cglab/human_chr21.gtf -N 5000 5000 -star_tmp /home/ghovhannisyan/TMP69 -t 7" ]


def getArgv(i):
    #print(testCases[i].split(' '))
    args = [] 
    for r in testCases[i].split(' '):
        if r.strip() != '':
            args.append(r)
            
    
    return args

#%% just test parser creatation
def testCreateArgumentParser():
    import crossmapper.crossmapper
    sys.argv =  getArgv(0)
    parser = crossmapper.crossmapper.createArgumentParser()
    parsedArgs = crossmapper.crossmapper.parseArgument(parser)
    #print(parsedArgs)
    return parsedArgs



def testRunCrossmap():
    sys.argv =  getArgv(1)
    #print(sys.argv)
    crossmapper.crossmapMain()

def testSimulation(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Simulation Test")
    crossmapper.mapping.concatGeneomes(parsedArgs)
    crossmapper.simulateReads.simulateData(parsedArgs)
    crossmapper.crossmapper.printOutputFileInfo(parsedArgs,'wgsim')
    logger.info("Simulation Test END")


def testMapping(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Mapping Test")
    crossmapper.mapping.mapping(parsedArgs)
    crossmapper.crossmapper.printOutputFileInfo(parsedArgs,'mapping')
    logger.info("Mapping Test END")    

def testCountingStep(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting Counting Test")

    res = crossmapper.countUtil.getReadCounters(parsedArgs)
    logger.info("Counting Test END")    
    return res
#testCreateArgumentParser()
#testRunCrossmap()

def renameChromosomes(parsedArgs):
    logger = crossmapper.helpers.getLogger()
    logger.info("Starting renameChromosomes Test")
    res = crossmapper.simulateReads.renameChromosomes(parsedArgs)
    logger.info("renameChromosomes Test END")    
    return res


#sys.argv =  getArgv(5)

args = testCreateArgumentParser()
renameChromosomes(args)
testSimulation(args)
testMapping(args)
res = testCountingStep(args)
#%%
# importlib.reload(crossmap)
#args = testCreateArgumentParser()
#res = testCountingStep(args)
# args.groupBarChart = False
# crossmap.reporting.createHTMLReport(res,args)
## write result




#%%
# importlib.reload(crossmap)
# res = testCountingStep(args)


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

