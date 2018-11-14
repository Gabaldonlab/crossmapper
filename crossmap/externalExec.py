## code for running external tools


import os
import sys
import subprocess

def execute(cmd,soft="def"):
    cmd_list = cmd.split(" ") 
    with open(f"{soft}_stdout.txt", "w") as outfile, open("{soft}_stderr.txt", "w") as errorfile:
        try:
            process = subprocess.run(cmd_list, shell=False, stdout = outfile, stderr = errorfile, check = False)
            print(process.returncode)
            if process.returncode != 0:
                #sys.exit("Command is wrong")
                print("Command is wrong")
                
        except FileNotFoundError as no_file:
            print(no_file, "NO SUCH FILE OR DIRECTORY")
        except Exception as ex:
            print(ex)
            #raise ex
            
            
        return


# =============================================================================
# def executeGffreadConvert(cmd_gffread_convert):
#     cmd_gffread_convert
#     
#     
# def executeGffreadExtract(cmd_gffread_extract):
#     cmd_gffread_extract
# 
# def executeWgsim(wgsim_cmd):
#     wgsim_cmd
# 
# def executeStarIndex(cmd_star_index):
#     cmd_star_index
#     
#     
# def executeBwaIndex(cmd_bwa_index):
#     cmd_bwa_index
#     
#     
# def executeStarMapping(cmd_star_mapping):
#     cmd_star_mapping
#     
# 
# def executeSamtoolsSort(cmd_samtools_sort):
#     cmd_samtools_sort    
#     
#     
#     
# def executeBwaMapping(cmd_bwa_mapping):
#     cmd_bwa_mapping
#     
# =============================================================================
    
    
    
    
