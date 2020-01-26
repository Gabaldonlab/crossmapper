#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:14:13 2019

@author: ahmed
"""

from enum import Enum
import re
from crossmapper.helpers import getBaseName, getLogger
import crossmapper.externalExec
from Bio import SeqIO
import math
import os 
## globals().update(dy['math'].__dict__)
## __import__("math")

#%%
class MapperType(Enum):
    DNA = 1
    RNA = 2
    BOTH = 3


class TemplateParseException(Exception):
    """Raised when can not parse command template"""
    pass

class ExecException(Exception):
    """Raised when can not execute the command"""
    pass
class NotImplException(Exception):
    pass
#%%
# {outputfile} : output bam/sam file , this whould be file name with full path {base_dir}/{mapper_name}_output/[internal pattern with readlen and layout e.g concat_50_PE  ].{output_type}
_var_outputfile = "outputfile"
_var_outputfile_prefix = "outputfile_prefix"
_var_output_type = "output_type"
_var_mapper_name = "mapper_name"
_var_outputfile_pattern = "outputfile_pattern"

# {gtf} : gtf file
_var_gtf = "gtf"
# {base_dir} : -o option of crossmaper
_var_base_dir = "base_dir"
# {out_dir} : output dir this will be the whole path {base_dir}/{mapper_name}_output    [-o option of crossmaper {base_dir} ] + outputdir option in the template {mapper_name} + prefix "_output"
_var_output_dir = "output_dir"
# {ref_dir} : reference index dir this will {base_dir}/{mapper_name}_index [-o option of crossmaper] + outputdir option in the template + postfix "_index"
_var_ref_dir  = "ref_dir"
# {ref_prefix} :: prefex name to the genome index
_var_ref_prefix = "ref_prefix"
# {ref_index} : {ref_dir}/{ref_prefix} 
_var_ref_index = "ref_index"
# {1} {2} : fastq files
_var_fastq1 = "1"
_var_fastq2 = "2"
# {ref_fasta} : reference fasta file
_var_ref_fasta = "ref_fasta"

# {layout} : PE or SE
_var_layout = "layout"
# {read_len} : read len 
_var_read_len  = "read_len"
_var_genome_len = "genome_len"

_var_n_threads = "n_threads"
_var_tmp_dir = "tmp_dir"



_var_qualifier_basename = "basename"


## [required]
_template_key_type = "type"
_template_key_dep = "dep"
_template_key_template = "template"
_template_key_index = "index"
_template_key_se = "se"
_template_key_pe = "pe"
_template_key_both = "both"

_template_key_mapper_name = _var_mapper_name
_template_key_sorted = "sorted"
_template_key_output_type = _var_output_type
_template_key_outputfile_pattern = _var_outputfile_pattern 


#%%
class Env():
    def __init__(self,args = None  ):
        self.args = args
        self.sharedEnv = {}
        self.qualifiers = {}
        
        self.qualifiers[_var_qualifier_basename] = qualifierBaseName

    def setEnvVar(self,varKey,varValue):
        if isinstance(varValue,str) :
            self.sharedEnv[varKey] = CVarToken(varValue)
        else:
            self.sharedEnv[varKey] = varValue
    def setQualifier(self,qualifierKey,qualifier):
        self.qualifiers[qualifierKey] = qualifier
    def getEnvVarValue(self,varKey):
        if varKey in self.sharedEnv:
            varValue = self.sharedEnv[varKey]
            if isinstance(varValue,CMDToken) :
                return varValue.eval(self)
            else:
                return varValue
        raise Exception("Can not evaluate unknown variable ({0}) used in the template. Please refer to the manual for list of available variables.".format(varKey))
        # raise Exception("Can not evaluate unknown variable ({0}) used in the template. Please refer to the manual for list of available variables.".format(varKey))

    def eval(self,varKey,qualifier = None):
        varValue = self.getEnvVarValue(varKey)
        if qualifier is None :
            return str(varValue)
        else:
            if qualifier in self.qualifiers :
                return str(self.qualifiers[qualifier](varKey,varValue,self))
            else :
                raise Exception("Unknown qualifier ({0}) used in the template. Please refer to the manual for list of available qualifier.".format(qualifier))


def checkMisssing(templateDict,key):
    if not key in templateDict:
        raise Exception('Mapper template missing required config key item ({0})'.format(key))
    return templateDict[key]

     
class CMDToken:
    
    def __init__(self,term):
        ## try to see if it work file
        self.term =  str(term)
        return

    def eval(self,env = None):
        return self.term
    def __repr__(self):
        return "\"{0}\"".format(self.term)

class LiteralToken(CMDToken):
    pass

TEMPLATE_SCANNER = re.compile(r'''
  (\s+) |                    # whitespace
  ([^{}\s<>\\*\\\(\\)]+) |   # Strings
  {((?:[^}\n\\]|\\.)*)} |    # variables
  >(.*) |                      # redirect
  (.)                        # an error!
''', re.DOTALL | re.VERBOSE) 

VAR_SCANNER = re.compile(r'''
  ([^{}\s<>\\*\\\(\\)]+) |   # Strings
  {((?:[^}\n\\]|\\.)*)} |    # variables
  (.)                        # an error!
''', re.DOTALL | re.VERBOSE)     

class VarToken(CMDToken):
    def __init__(self,term):
        term = str(term)
        super().__init__(term)
        if ":" in term:
            tokens = term.split(":")
            self.varName = tokens[0]
            self.qualifier = tokens[1]
        else:
            self.varName = term
            self.qualifier = None
    
    
    def __repr__(self):
        if self.qualifier:
            return "\"{0}:{1}\"".format(self.varName,self.qualifier)
        return "\"{0}\"".format(self.varName)
    
    def eval(self,env):
        if isinstance(env, Env) :
            return env.eval(self.varName,self.qualifier)
        else:
            raise Exception("Invalid method Args. Unknown Env Object in eval")

class CVarToken(CMDToken):
    def __init__(self,term):
        super().__init__(term)
        self.varTemplate = str(term)
        
        self.varTemplateTokens = []
        
        for match in re.finditer(VAR_SCANNER, self.varTemplate):
            literalToken,variableToken,errors =  match.groups()
            # print([literalToken,variableToken,errors])
            if errors :
                s,e = match.span()
                errMsg = "Unexpecet token or char highlighted between <<?>> in the following VAR : {0}<<{1}>>{2}".format(self.varTemplate[:s],self.varTemplate[s:e],self.varTemplate[e:])
                raise Exception(errMsg)
            if variableToken:
                self.varTemplateTokens.append(VarToken(variableToken))
            elif literalToken:
                self.varTemplateTokens.append(LiteralToken(literalToken))

    
    def __repr__(self):
        return "\"{0}\"".format(self.varTemplate)
    
    def eval(self,env):
        if isinstance(env, Env) :
            varValue = ""
            for tkn in self.varTemplateTokens:
                varValue = varValue + tkn.eval(env)
            return varValue
        else:
            raise Exception("Invalid method Args. Unknown Env object in eval")
   


def qualifierBaseName(varName,varValue,env = None):
    # varValue = cmdTemplate.env[inVar]
    return os.path.splitext(os.path.basename(varValue))[0]


## main class for a command template
class CMDTemplate:
    
    def __init__(self,template, env = None , cmdName = None):
        
        self.template = template
        
        self.cmdName = cmdName
        self.env = env
       
        self.redirectIndex = -1
        self.stdoutRedir = None
        self.parse()
        
        
    def parse(self):
        ## first ckeck if we had > operator inside the command
        
#        if ">" in self.template :
#            ## for now assume only > one operator stdout forward to file 
#            tokens = self.template.split(">")
#            self.templateBase = tokens[0].strip()
#            self.stdoutRedir = tokens[1].strip()
#        else:
#            self.templateBase = self.templateBase = tokens[0].strip() 
        
        self.templateTokens = []
        # print(self.template)
        for match in re.finditer(TEMPLATE_SCANNER, self.template):
            whitespaces,literalToken,variableToken,redir,errors =  match.groups()
            # print([whitespaces,literalToken,variableToken,redir,errors])
            if errors :
                s,e = match.span()
                errMsg = "Unexpecet token or char highlighted between <<?>> in the following line :\n{0}<<{1}>>{2}".format(self.template[:s],self.template[s:e],self.template[e:])
                raise Exception(errMsg)
            if variableToken:
                self.templateTokens.append(VarToken(variableToken))
            elif redir :
                self.redirectIndex = len(self.templateTokens)
                self.stdoutRedir = CVarToken(redir.strip())
            elif whitespaces:
                self.templateTokens.append(LiteralToken(whitespaces))
            elif literalToken:
                self.templateTokens.append(LiteralToken(literalToken))

    def setEnv(self,env):
        self.env = env
    
    def execCmd(self, cmdDir = None):
        if cmdDir is None :
            # get the cmdDir from the env 
            # this should be 
            cmdDir =  self.env.getEnvVarValue (_var_output_dir)
        if self.cmdName is None:
            ## the first token should be the name
            cmdName = self.templateTokens[0].eval(self.env)
        else:
            cmdName = self.cmdName
            
        cmdTemplateStr = ""
        for term in self.templateTokens:
            cmdTemplateStr += term.eval(self.env)
            
        if self.stdoutRedir is None :
            res = crossmapper.externalExec.execute(cmdTemplateStr, f"{cmdName}" , outDir = f"{cmdDir}", overwrite = False)
            if not res.resCheck(clean = False):
                raise ExecException(f"Failed to excute command {cmdTemplateStr},\nSee log file {res.stdErrFile}")
        else:
            redirectFile = self.stdoutRedir.eval(self.env)
            res = crossmapper.externalExec.execute(cmdTemplateStr, f"{cmdName}", redirectFile,None,  outDir = f"{cmdDir}")
            if not res.resCheck():
                raise ExecException(f"Failed to excute command {cmdTemplateStr},\nSee log file {res.stdErrFile}")
                

class CMDTemplateList:
    def __init__(self,templateList, env = None , cmdName = None):
        
        if isinstance(templateList,str):
            self.templateList = [templateList]
        else:
            self.templateList = templateList
        
        self.cmdName = cmdName
        self.env = env
        self.cmdTempaltes = []
        
        for template in self.templateList:
            self.cmdTempaltes.append(CMDTemplate(template,env,cmdName))
        
        
    def setEnv(self,env):
        self.env = env
    def execCmdList(self, cmdDir = None):
        for cmdTemplate in self.cmdTempaltes:
            cmdTemplate.setEnv(self.env)
            cmdTemplate.execCmd(cmdDir)

#%% #################################################
class MapperBase:
    '''base class for all Mapper'''
    def __init__(self , parsedArgs, mapperName = "mapper" ):
        self.parsedArgs = parsedArgs
        self.logger = getLogger()
        self.mapperName = mapperName
        self.sorted = False
        self.indexFolder = f"{parsedArgs.out_dir}/{mapperName}_index"
        self.mappingDir = f"{parsedArgs.out_dir}/{mapperName}_output"
    def run(self):
        
        ## TODO: :what to add here
        
        ## TODO: check if every thing is setup
        
        
        self.runBuildIndex()
        self.runMapping()
        pass
    def runBuildIndex(self):
        # calcualte concat.fasta genome size
        '''
        self.indexFolder must be set before calling this method, or default value will be used
        '''
        self.genomeLen = 0
        for rec in SeqIO.parse(self.parsedArgs.genomeConcatFasta, 'fasta'):
            self.genomeLen+=len(rec.seq)
            
        if self.genomeLen > 3000000000:
            self.logger.warning("Concatenated genome size is larged than 3 Gb!." )
        
        
        #SA_index_size = min(14, round(math.log(genome_len,2)/2) - 1)
        #logger.debug("genomeSAindexNbases = %s"%(SA_index_size))
        self.logger.info(f"Starting genome indexing for {self.mapperName}.")
        
        
        if os.path.isdir(f"{self.indexFolder}") == True:
            self.logger.info(f"{self.indexFolder} directory exists. Generating index files.")
        else:
            self.logger.info(f"Creating {self.indexFolder} directory. Writing index files to {self.indexFolder}.")
            os.makedirs(f"{self.indexFolder}")
        
        self.onRunBuildIndex(self.genomeLen,self.indexFolder)
        
        self.logger.info(f"Genome index for {self.mapperName} is generated.")

        pass
    
    def runMapping(self):
        '''
        Do not Override, if so make sure to call super().runMapping() after subclass setup
        self.mappingDir must be set before calling this method, or defaul value will be used
        Method also use
        self.mapperName
        '''

        self.parsedArgs.mappingDir =  self.mappingDir
        
       
        if os.path.isdir(f"{self.mappingDir}") == False:
            self.logger.info(f"Creating {self.mappingDir} directory.")
            os.makedirs(f"{self.mappingDir}")
        self.logger.info(f"Starting {self.mapperName} Mapper.")

        
        
        for rlen in self.parsedArgs.input_rlen:
            self.parsedArgs.mappingOutputFiles[rlen] = {}
            
            if self.parsedArgs.read_layout != "PE" : ## i.e SE or both
                # Single end case
                layout = "SE"
                tmpFile = self.doMap(self.parsedArgs.simulationOutputFiles[rlen],rlen,layout)
                finalBamFile = f"{self.mappingDir}/concat_{rlen}_{layout}_sorted.bam"
                if self.sorted :
                    finalBamFile = tmpFile
                else:
                    self.sortBam(tmpFile,finalBamFile,self.mappingDir)
                self.indexBam(finalBamFile,self.mappingDir)
                self.parsedArgs.mappingOutputFiles[rlen][layout] = finalBamFile
                ####
            if self.parsedArgs.read_layout != "SE": ## i.e. PE or both
                layout = "PE"
                tmpFile = self.doMap(self.parsedArgs.simulationOutputFiles[rlen],rlen,layout)
                finalBamFile = f"{self.mappingDir}/concat_{rlen}_{layout}_sorted.bam"
                if self.sorted :
                    finalBamFile = tmpFile
                else:
                    self.sortBam(tmpFile,finalBamFile,self.mappingDir)
                self.indexBam(finalBamFile,self.mappingDir)
                self.parsedArgs.mappingOutputFiles[rlen][layout] = finalBamFile

        pass

##########################  Must be implemented in subclasess
    def doMap(self,fastqFiles, rlen, read_layout):
        raise NotImplException("Sub Class must implement doMap")
    def onRunBuildIndex(self,genomeLen,indexFolder):
        raise NotImplException("Sub Class must implement onRunBuildIndex")
##########################  Must be implemented in subclasess

    def sortBam(self,inFile,outFile,runDir,deleteInFile = True):
        res = crossmapper.externalExec.execute(f"samtools sort -@{self.parsedArgs.threads} -o {outFile} {inFile}",
                                        "samtools" ,
                                        outDir = f"{runDir}")
        if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
            raise ExecException("samtools sort Execution Failed.")
        
        self.logger.info("Sorting is finished. " +f"Final bam file writen to {outFile}")

        if deleteInFile:
            try :
                self.logger.debug(f"Deleteing tmp sam file {inFile}")
                os.remove(inFile)
            except :
                self.logger.warning(f"Can not delete temporary bam file {inFile}")
            #sys.exit("Execution Failed")

    def indexBam(self,inBamFile, runDir):
        cmd_samtools_index = f"samtools index {inBamFile}"
        # print(cmd_samtools_index)
        res = crossmapper.externalExec.execute(cmd_samtools_index,"samtools_index" , outDir = f"{runDir}")
    
        if not res.resCheck(stdoutRemove=True,stdErrRemove=True):
            raise ExecException("samtools index Execution Failed")
            # sys.exit("Execution Failed")
       
    
class STARMapper(MapperBase):
    def __init__(self,parsedArgs):
        super().__init__(parsedArgs,"STAR")
    def onRunBuildIndex(self,genomeLen,  indexFolder ):
        if genomeLen > 3000000000:
            self.logger.warning("STAR mapping will require more than 30 GB of RAM." )
    
    
        SA_index_size = min(14, round(math.log(genomeLen,2)/2) - 1)
        self.logger.debug("genomeSAindexNbases = %s"%(SA_index_size))

        cmd_star_index = "STAR --runMode genomeGenerate " \
        f"--runThreadN {self.parsedArgs.threads} " \
        f"--genomeDir {indexFolder} " \
        f"--genomeFastaFiles {self.parsedArgs.genomeConcatFasta} " \
        f"--genomeSAindexNbases {SA_index_size}"
    
        res = crossmapper.externalExec.execute(cmd_star_index,"STAR_index" , outDir = f"{self.parsedArgs.out_dir}" )
        if not res.resCheck( stdoutRemove = True ):
            raise ExecException("STAR Index Execution Failed")
            #sys.exit("Execution Failed")

    def doMap(self,fastqFiles, rlen, read_layout):
        #reads =   " ".join(fastqFiles).strip()
        if read_layout == "PE" :
            reads =   " ".join(fastqFiles).strip()
        else:
            reads =   fastqFiles[0]
        overhang = rlen - 1
        star_dir = self.mappingDir
        if self.parsedArgs.bacterial_mode is True:
            intron_len_max=1
        else:
            intron_len_max=0
                
        cmd_star_mapping = "STAR " \
            f"--runThreadN {self.parsedArgs.threads} " \
            f"--genomeDir {self.indexFolder} " \
            f"--sjdbGTFfile {self.parsedArgs.out_dir}/concat.gtf " \
            f"--sjdbOverhang {overhang} " \
            f"--readFilesIn {reads} " \
            "--readFilesCommand cat --outSAMtype BAM Unsorted " \
            f"--outFileNamePrefix {star_dir}/concat_{rlen}_{read_layout}_ " \
            f"--outFilterMismatchNmax {self.parsedArgs.outFilterMismatchNmax} " \
            f"--outFilterMultimapNmax 10000 " \
            f"--outFilterMismatchNoverReadLmax {self.parsedArgs.outFilterMismatchNoverReadLmax} " \
            f"--alignIntronMax {intron_len_max} " \
            f"--outTmpDir {self.parsedArgs.star_temp_dir}"
        
        # print(cmd_star_mapping)
        res = crossmapper.externalExec.execute(cmd_star_mapping,"STAR_mapping" , outDir = f"{star_dir}", overwrite = False)
        if not res.resCheck(clean = False):
            raise ExecException("STAR Execution Failed.")
            # sys.exit("Execution Failed")
        tmpBamFile = f"{star_dir}/concat_{rlen}_{read_layout}_Aligned.out.bam"
        return tmpBamFile

class BWAMapper(MapperBase):
    def __init__(self,parsedArgs):
        super().__init__(parsedArgs,"BWA")
    def onRunBuildIndex(self,genomeLen,  indexFolder ):
        algo = ""
        if genomeLen > 3000000000:
            self.logger.info("Using bwtsw algorithm for index generation." )
            algo = "-a bwtsw"
        
        
        cmd_bwa_index = "bwa index " \
        f"-p {indexFolder}/concat_BWA " \
        f"{algo} " \
        f"{self.parsedArgs.genomeConcatFasta}"

        # print(cmd_bwa_index)
        res = crossmapper.externalExec.execute(cmd_bwa_index,"BWA_index", outDir = f"{self.parsedArgs.out_dir}")
        if not res.resCheck( stdoutRemove = True ):
            raise ExecException("bwa index Execution Failed.")
            # sys.exit("Execution Failed")
        
    pass

    def doMap(self,fastqFiles,rlen,read_layout):
        if read_layout == "PE" :
            reads =   " ".join(fastqFiles).strip()
        else:
            reads =   fastqFiles[0]
        bwa_dir = self.mappingDir 
        
        tmpSamFile = f"{bwa_dir}/concat_{rlen}_{read_layout}.sam"
        
        cmd_bwa_mapping = f"bwa mem -a -t {self.parsedArgs.threads} -A {self.parsedArgs.match_score} -B {self.parsedArgs.mismatch_penalty} {self.parsedArgs.out_dir}/BWA_index/concat_BWA {reads}" 
        #f"samtools sort @{parsedArgs.threads} -o concat_{rlen}_{read_layout}_sorted.bam -
        res = crossmapper.externalExec.execute(cmd_bwa_mapping,"BWA_mapping" , 
                                     tmpSamFile,
                                     None,
                                     outDir = f"{bwa_dir}")
        if not res.resCheck():
            raise ExecException("bwa Execution Failed.")
            #sys.exit("Execution Failed")
        # finalBamFile = f"{bwa_dir}/concat_{rlen}_{read_layout}_sorted.bam"
        return tmpSamFile

class TemplateMapper(MapperBase):
    'General class for mapper tools'
    def __init__(self , configTemplate,parsedArgs ):
        '''
        template is dict parsed from a yaml config file.
        '''
        super().__init__(parsedArgs)
        # self.parsedArgs = parsedArgs

        self.configTemplate = configTemplate
        ## setting default values
        self.sorted = False
        self.type = MapperType.BOTH
        self.outputType = "sam"
        self.mapperName = "mapper"
        self.output = "{out}." + self.outputType
        self.dep = None
        
        ## optional
        ## check 
        if _template_key_type in configTemplate :
            ## type could be RNA/DNA or Both
            typeOptions = ['DNA','RNA','Both']
            typeValue = configTemplate[_template_key_type]
            if typeValue in typeOptions:
                self.type = MapperType(typeOptions.index(typeValue) + 1)
            else:
                raise Exception('Unknown type option value ({0}) in the mapper config template'.format(typeValue))
        if _template_key_output_type in configTemplate :
            outputTypeOptions = ['sam','bam']
            outputTypeValue = configTemplate[_template_key_output_type]
            if outputTypeValue in outputTypeOptions:
                self.outputType = outputTypeValue
            else:
                raise Exception('Unknown output_type option value ({0}) in the mapper config template. Valid options are sam/bam.'.format(outputTypeValue))
        
        if _template_key_sorted in configTemplate :
            sortedValue = configTemplate[_template_key_sorted]
            if type(sortedValue) == bool :
                self.sorted = sortedValue
            else:
                raise Exception('Unknown sorted option value ({0}) in the mapper config template. Valid options are yes/no.'.format(sortedValue))
        
        if _template_key_mapper_name in configTemplate :
            self.mapperName = configTemplate[_template_key_mapper_name]
            if re.search(r'[^A-Za-z0-9_\-\\]',self.mapperName): 
                raise Exception('mapper_name config ({0}) contains invalid chars'.format(self.mapperName))
        
        if _template_key_outputfile_pattern in configTemplate :
            self.outputFilePattern = configTemplate[_template_key_outputfile_pattern]
        else:
            self.outputFilePattern = f"{{{_var_outputfile_prefix}}}.{{{_var_output_type}}}"

        if _template_key_dep in configTemplate:
            self.dep = configTemplate[_template_key_dep]
            tokens = self.dep.split("/")
            if len(tokens) > 1 :
                self.channelName = tokens[0]
                self.packageName = tokens[1]
            else:
                self.channelName = None
                self.packageName = tokens[1]

        ## required
        self.template = checkMisssing(configTemplate,_template_key_template)
        ## template should have index
        self.indexTemplate = checkMisssing(self.template,_template_key_index)
        ## check if has both
        if _template_key_both in self.template:
            self.seTemplate = self.template[_template_key_both]
            self.peTemplate = self.template[_template_key_both]
        else:
            self.seTemplate = checkMisssing(self.template,_template_key_se)
            self.peTemplate = checkMisssing(self.template,_template_key_pe)


    def checkDep(self):
        if not self.dep is None:
            try :
                checkCmd = f'conda list | grep "{self.packageName}" -w | wc -l'
                installCmd = f'conda install {self.packageName} -y'
                if not self.channelName is None :
                    installCmd = f"{installCmd} -c {self.channelName}"
                import subprocess
                checkProcess = subprocess.run(checkCmd,shell=True, stdout = subprocess.PIPE)
                checkRes = checkProcess.stdout
                if isinstance(checkRes,bytes):
                    checkRes = checkRes.decode("UTF-8")
                if isinstance(checkRes,str):
                    checkRes = checkRes.strip()
                    if not checkRes == "1" :
                        self.logger.info(f"conda package {self.packageName} is not installed, Trying to install it")
                        ## package does not exist, install
                        installProcess = subprocess.run(installCmd,shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                        if installProcess.returncode == 0:
                            self.logger.info(f"Finished installing conda package {self.packageName}.")
                        else:
                            raise Exception("Conda install Fail.\n" + installProcess.stderr.decode())
            except Exception as ex:
                raise Exception(f"Can not conda install mapper package {self.packageName}, Please try to install it manully before running crossmapper. " + str(ex))
            
    def run(self):
        
        # setup shared env variables
        self.sharedEnv = Env(self.parsedArgs)
        
        
        self.sharedEnv.setEnvVar(_var_tmp_dir,self.parsedArgs.star_temp_dir)
        self.sharedEnv.setEnvVar(_var_n_threads,self.parsedArgs.threads)
        
        # self.sharedEnvVar[_var_base_dir] = self.parsedArgs.out_dir
        self.sharedEnv.setEnvVar(_var_base_dir,self.parsedArgs.out_dir)
        
        self.sharedEnv.setEnvVar(_var_mapper_name,self.mapperName)
        
        self.sharedEnv.setEnvVar(_var_output_dir, "{base_dir}/{mapper_name}_output" )
        
        
        ## ref
        self.sharedEnv.setEnvVar(_var_gtf , self.parsedArgs.out_dir +"/concat.gtf")
        self.sharedEnv.setEnvVar(_var_ref_dir,f"{{{_var_base_dir}}}/{{{_var_mapper_name}}}_index")
        self.sharedEnv.setEnvVar(_var_ref_prefix, "concat_index")
        self.sharedEnv.setEnvVar(_var_ref_index,   f"{{{_var_ref_dir}}}/{{{_var_ref_prefix}}}")
        self.sharedEnv.setEnvVar(_var_ref_fasta, self.parsedArgs.genomeConcatFasta)
        
        self.sharedEnv.setEnvVar(_var_outputfile, f"{{{_var_output_dir }}}/{{{_var_outputfile_pattern}}}")
        self.sharedEnv.setEnvVar(_var_outputfile_pattern, self.outputFilePattern)
        self.sharedEnv.setEnvVar(_var_output_type,self.outputType)
        
        
        
        
        
        ## setup needed variables values
        
        
        
        self.indexFolder = self.sharedEnv.getEnvVarValue(_var_ref_dir)
        self.mappingDir = self.sharedEnv.getEnvVarValue(_var_output_dir)
        
        
        self.seCmdTemplate = CMDTemplateList(self.seTemplate)
        self.peCmdTemplate = CMDTemplateList(self.peTemplate)

        self.seCmdTemplate.setEnv(self.sharedEnv)
        self.peCmdTemplate.setEnv(self.sharedEnv)
        
        ## first run build index template
        super().run()
        
        ## for all input files start mapping cmd
        pass
    
    def onRunBuildIndex(self,genomeLen,indexFolder):
        
        indexCmdTemplate = CMDTemplateList(self.indexTemplate,env = self.sharedEnv,  cmdName = self.indexFolder )        
        indexCmdTemplate.execCmdList(cmdDir = self.indexFolder)
        
                
    def doMap(self,fastqFiles,rlen,read_layout):
        
        
        
        self.sharedEnv.setEnvVar(_var_read_len,rlen)
        self.sharedEnv.setEnvVar(_var_fastq1,fastqFiles[0])
        self.sharedEnv.setEnvVar(_var_fastq2,"")
        if read_layout == 'SE':
            # Single end case
            self.sharedEnv.setEnvVar(_var_outputfile_prefix,f"concat_{rlen}_SE")
            self.sharedEnv.setEnvVar(_var_layout,"SE")
            # se_mapping = fastqFiles[0]
            # print(se_mapping,rlen,"SE")
            self.seCmdTemplate.execCmdList()
        if read_layout == 'PE' :
            self.sharedEnv.setEnvVar(_var_outputfile_prefix,f"concat_{rlen}_PE")
            self.sharedEnv.setEnvVar(_var_fastq2,fastqFiles[1])
            self.sharedEnv.setEnvVar(_var_layout,"PE")
            #pe_mapping = " ".join(fastqFiles)
            # print(pe_mapping,rlen,"PE")
            self.peCmdTemplate.execCmdList()
        return self.sharedEnv.getEnvVarValue(_var_outputfile)

