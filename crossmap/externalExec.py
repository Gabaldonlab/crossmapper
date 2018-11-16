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
        except PermissionError as perm_denied:
            print(perm_denied, "PERMISSION DENIED")
        except Exception as ex:
            print(ex)
            #raise ex
            
            
        return
