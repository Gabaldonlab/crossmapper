import pysam
import sys
import os
from crossmap.helpers import getLogger


 
#%%

def getSequencesPerOrganism(inFastaFile):
    '''
        get/read sequence names from fastq file for one Genome
    '''
    inFile = open(inFastaFile,"r")
    seqsDir = []
    for line in inFile:
        if line.startswith(">"):
            line = line.replace(">","").strip()
            line = line.split(" ")
            seqsDir.append(line[0])
    return seqsDir
def getSequencesPerOrganisms(listFiles):
    '''
        read fasta files and return all sequences names
    '''
    seqs = []
    for oneFile in listFiles:
       seqs.append( getSequencesPerOrganism(oneFile))
    return seqs


############## Section to read gtf file and get/map transcript names to sequences name in the genome 
def getAtts(attToken):
    ## for gtf sep is ;
    attsDict= {}
    tokens = attToken.strip().split(";")
    # print(attToken)
    for oneAtt in tokens :
        ## sep is space for now
        if oneAtt == "":  # if empty just continue
            continue
        
        ## TODO check if quated then whatch/scan for inner spaces 
        keyValue = oneAtt.strip().split(' ')
        
        
        attsDict[keyValue[0]] = keyValue[1].replace("\"","").strip()
    return attsDict


def parseGTFrecord(record, selectedFeatures = None):
    tokens = record.split("\t")
    ## placeholder :: TODO parse everything here if needed
    ## TODO :: check GTF record if it is OK
    ## warn usr if somethng is wrong with this line
    recordDict = {}
    
    recordDict['seqName'] = tokens[0]
    
    
    
    attsToken = tokens[8]
    recordDict.update(getAtts(attsToken))
    if selectedFeatures != None :
        return dict((k, recordDict[k]) for k in selectedFeatures)
    return recordDict
    

def mapTranscriptToSequence(inGTFFiles , seqsIndex = None):
    transcriptMap = {}
    
    for oneGTFFile in inGTFFiles:
        
        inFile = open(oneGTFFile,'r')
        for line in inFile:
            line = line.strip()
            recordDict = parseGTFrecord(line, ['seqName','transcript_id'])
            if seqsIndex != None :
                transcriptMap[ recordDict['transcript_id'] ] = seqsIndex[recordDict['seqName']]
            else:
                transcriptMap[ recordDict['transcript_id'] ] = recordDict['seqName']
            # seqName,transcript_id = 
    return transcriptMap
##%%%%%%%%%%%%
def mapTranscriptToSequenceGFF(inGFFFiles):
    transcriptMap = {}
    
    for oneGTFFile in inGFFFiles:
        
        inFile = open(oneGTFFile,'r')
        for line in inFile:
            line = line.strip()
            recordDict = parseGTFrecord(line, ['seqName','ID'])
            transcriptMap[ recordDict['ID'] ] = recordDict['seqName']
            # seqName,transcript_id = 
    return transcriptMap


#transcriptMapGFF = mapTranscriptToSequenceGFF([sp1InputGFF])

#%%
# mock, TODO 
#def getSequencesPerOrganism(speciesIds):
#    seqs = {}
#    seqs[speciesIds["C_parapsilosis_CDC317"]] = []
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig005504_C_parapsilosis_CDC317") 
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig005569_C_parapsilosis_CDC317")
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig005806_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig005807_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig005809_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig006110_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig006139_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("Contig006372_C_parapsilosis_CDC317")   
#    seqs[speciesIds["C_parapsilosis_CDC317"]].append("mito_C_parapsilosis_CDC317")
#    
#    seqs[speciesIds["C_albicans_SC5314"]] = []
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr1A_C_albicans_SC5314")
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr2A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr3A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr4A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr5A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr6A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chr7A_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chrM_C_albicans_SC5314")  
#    seqs[speciesIds["C_albicans_SC5314"]].append("Ca22chrRA_C_albicans_SC5314") 
#    return seqs

def sequenceToOrganism(allSeqs):
    seqToOrg = {}
    for orgIndex in range(0,len(allSeqs)):
#    for org,seqs in allSeqs.items():
        for seq in allSeqs[orgIndex] :
            seqToOrg[seq] = orgIndex
    return seqToOrg

def createSequenceIndex(allSeqs):
    seqIndex = {}
    curI = 0
    for seqList in allSeqs:
        for seq in seqList:
            seqIndex[seq] = curI
            curI +=1
    return seqIndex

#%%
    
#%%
def parseWgsimQNAME(qName,isPairEnd=True):
    ## xyz_start_end_int:int:int_int:int:int_hexa
    tokens = qName.split('_')
    memberPairIndex = tokens[-1] # last token
    rightPairInfo = ''
    rightPairInfo = ''
    if isPairEnd:
        rightPairInfo = tokens[-2]
        leftPairInfo  = tokens[-3]
    else:
        ## TODO :: I can not see any difference in how wgsim handle ids 
        rightPairInfo = tokens[-2]
        leftPairInfo  = tokens[-3]
    end = int(tokens[-4])
    start = int(tokens[-5])
    orgSeqName = "_".join(tokens[:-5])
    readId =  qName # "_".join(tokens[:-3])
    # print("%s \t%s\t%d\t%d\t%s\t%s" % (qName, orgSeqName, start, end , leftPairInfo,rightPairInfo ))

    return {    'readId': readId,
                 'orgSeqName' : orgSeqName,
                 'start' : start,
                 'end' : end,
                 'leftPairInfo': leftPairInfo,
                 'rightPairInfo': rightPairInfo,
                 'memberPairIndex': memberPairIndex
            }
#%%
class ReadId():
    def __init__(self,readId,orgSeqName,start,end,**kwargs):
        self.readId = readId
        self.orgSeqName = orgSeqName
        self.start = int(start)
        self.end = int(end)
    def __repr__(self):
        return "%s\t%s\t%d\t%d" % (self.readId, self.orgSeqName, self.start, self.end  )
class WGSIMRead(ReadId):
    def __init__(self,*args,**kwargs):
        '''
        readId,orgSeqName,start,end,leftPairInfo,rightPairInfo,memberPair
        or QNAME
        '''
        
        self.leftPairInfo = ''
        self.rightPairInfo = ''
        self.memberPair = ''
        
        ## if given QNAME instead of args do parse it here
        if 'qName' in kwargs:
            kwargs = parseWgsimQNAME(kwargs['qName'])
        
        super(WGSIMRead,self).__init__(*args,**kwargs)
        if 'leftPairInfo' in kwargs:
            self.leftPairInfo = kwargs['leftPairInfo']
        if 'rightPairInfo' in kwargs:
            self.rightPairInfo = kwargs['rightPairInfo']
        if 'memberPair' in kwargs:
            self.memberPair = kwargs['memberPair']
    def __repr__(self):
        return "%s\t%s\t%s\t%s" % (self.readId, self.leftPairInfo, self.rightPairInfo, self.memberPair  )
          
    
#%% counters 
    
class Counter():
    ## Counters Unique Mapping
    ## #species 
    ## 0          [ #CorrectMapping, #unCorrectMappingCorrectSp, #unCorrectSp ]
    ## 1          [ #CorrectMapping, #unCorrectMappingCorrectSp, #unCorrectSp ]
    def __init__(self,speciesIds):
        nSpecies = len(speciesIds)
        self.speciesIds = speciesIds
        self.unique         = [[0] * 3 for i in range(nSpecies)]
        self.unmapped       = [0]  * nSpecies
        self.total       = [0]  * nSpecies
        self.mPrimary       = [[0] * 3 for i in range(nSpecies)]
        self.mSecondary     = [[0] * 3 for i in range(nSpecies)]
        ## Counter is  [#readsMultiMaptoOnlyOneSp,#readsMultiMaptoBothSp, #readsMultiMaptoOneWrongSp]
        self.multiReads     = [[0] * 3 for i in range(nSpecies)]
        self.crossSpMulti = [[0] * nSpecies for i in range(nSpecies)]
        self.crossSpUnique = [[0] * nSpecies for i in range(nSpecies)]
        self.crossSpTotal = [[0] * nSpecies for i in range(nSpecies)]
    def summary(self , outFile = sys.stdout):
        
        lpad = 25
        
        
        
        print("Summary Report \n" + "".join( "=" * 15 )  , file=outFile)
        print(      "\t{0}:\t{1}"
                      .format("Totol Reads".ljust(lpad) , sum(self.total) + sum(self.unmapped)),
                      file=outFile )
    
        totalUniqueMapped = 0
        for uniqueCounter in self.unique:
            totalUniqueMapped += sum(uniqueCounter)
        print(      "\t{0}:\t{1}"
                      .format("Unique Mapped Reads".ljust(lpad) , totalUniqueMapped),    
                      file=outFile)
        
        totalMultiMapped = 0
        for multiCounter in self.multiReads:
            totalMultiMapped += sum(multiCounter)
        print(      "\t{0}:\t{1}"
                      .format("Multi Mapped Reads".ljust(lpad) , totalMultiMapped),
                      file=outFile)
        print(      "\t{0}:\t{1}"
                      .format("UnMapped Mapped Reads".ljust(lpad) , sum(self.unmapped)) ,   
                  file=outFile)
        
        print("\t-----------------\n",file=outFile)
        #print("Cross Species Stats \n")
        #UniqeMapped
        print( "{0}\t{1}\t\t|\t\t\t{2}\t".format(" ".ljust(4),  "Correct","unCorrect"   ) ,file=outFile     )
        print( "{0}\t{1}\t\t|\t\t\t{2}\t".format(" ".ljust(4),  "-------","-------"   )   ,file=outFile  )
        ## Multi_org+other mapped to source plus map to some other sp.
        ## Multi_onlyOther map to any of other sp and does not map to source sp
        print( "{0}\t{1}\t{2}\t|\t{3}\t{4}\t{5}\t{6}".format(" ".ljust(4), "Unique" ,"Multi" ,   "Unique", "Unique_Cross", "Multi_org+other" ,  "Multi_onlyOther" )   ,file=outFile)
        for sp in self.speciesIds:
            spID = self.speciesIds[sp]
            #print(spID)
            print( "{0}\t{1}\t{2}\t|\t{3}\t{4}\t{5}\t{6}".format(sp.ljust(4),
                  self.unique[spID][0] ,
                  self.multiReads[spID][0] ,   
                  
                  self.unique[spID][1],
                  self.unique[spID][2],
                  self.multiReads[spID][1],  
                  self.multiReads[spID][2]),
                file=outFile)
        
            
        print("\t-----------------\n",file=outFile)
        print("Unique_Cross\n\t",end="",file=outFile)
        for sourceSp in self.speciesIds:
            print(sourceSp,end="\t",file=outFile)
        print("\n",file=outFile)
        for sourceSp in self.speciesIds:
            sourceSpID = self.speciesIds[sourceSp]
            print(sourceSp,end="\t",file=outFile)
            for targetSp in self.speciesIds:
                targetSpID = self.speciesIds[targetSp]
                print(self.crossSpUnique[sourceSpID][targetSpID],end="\t",file=outFile)
            print("\n",file=outFile)
        
        
        print("\t-----------------\n",file=outFile)
        print("Multi_Cross\n\t",end="",file=outFile)
        for sourceSp in self.speciesIds:
            print(sourceSp,end="\t",file=outFile)
        print("\n",file=outFile)
        for sourceSp in self.speciesIds:
            sourceSpID = self.speciesIds[sourceSp]
            print(sourceSp,end="\t",file=outFile)
            for targetSp in self.speciesIds:
                targetSpID = self.speciesIds[targetSp]
                print(self.crossSpMulti[sourceSpID][targetSpID],end="\t",file=outFile)
            print("\n",file=outFile)
        
        print("Total_Cross\n\t",end="",file=outFile)
        for sourceSp in self.speciesIds:
            print(sourceSp,end="\t",file=outFile)
        print("\n",file=outFile)
        for sourceSp in self.speciesIds:
            sourceSpID = self.speciesIds[sourceSp]
            print(sourceSp,end="\t",file=outFile)
            for targetSp in self.speciesIds:
                targetSpID = self.speciesIds[targetSp]
                print(self.crossSpTotal[sourceSpID][targetSpID],end="\t",file=outFile)
            print("\n",file=outFile)
                

#%% Read RNA
def checkNHTag(bamFile):
    
    for record in bamFile.fetch(until_eof=True):
        if record.has_tag('NH'):
            #
            return None
        else:
            break
    #print("No")
    bamFile.reset()
    allReads = {}
    for record in bamFile.fetch(until_eof=True):
        if record.is_unmapped :
            allReads[record.qname] = 0
            continue

        
        if record.is_proper_pair :
            if not record.qname in allReads :
                allReads[record.qname] = 0
            if record.is_paired and record.is_read1 :
                allReads[record.qname]+=1
        else:
            pairName =  record.qname+"_1" if record.is_read1 else   record.qname+"_2"
            if not pairName in allReads :
                allReads[pairName] = 0
            allReads[pairName]+=1
    bamFile.reset()
    return allReads
def countReads(bamFile, speciesIds , seqsIndex , seqToOrg , transcriptMap = None , nhTag = None):
    nSpecies = len(speciesIds)
    #uniqueCounters = [[0] * 3 for i in range(nSpecies)]
    #unmappedCounters = [0] * len(nSpecies)
    #mPrimaryCounters    =  [[0] * 3 for i in range(nSpecies)]
    #mSecondaryCounters  =  [[0] * 3 for i in range(nSpecies)]
    allCounter = Counter(speciesIds)
    count  = 0
    dict_list=[]
    reads = {}
    multiReadsClass = {}
    multiReadsMappedTpSp = {}
    multiReadsOrgSp = {}
    bamRecords = {}
    
    records = []
    
    for record in bamFile.fetch(until_eof=True):
        count+=1
        #if count == 100 :
        #    break
        #print(record.qname)
        read = WGSIMRead(qName=record.qname)
        #print(read)
        if transcriptMap != None:
            read.orgSeqName = seqsIndex[transcriptMap[read.orgSeqName]]
        
        ## count unmapped reads
        orgReadSpId = seqToOrg[read.orgSeqName ]
        if record.is_unmapped:
            allCounter.unmapped[orgReadSpId] +=1
            continue
        allCounter.total[orgReadSpId] +=1
        # read = WGSIMRead(qName=record.qname)
        mappingSpId = seqToOrg[record.reference_name]
        
        #print(mappingSpId,orgReadSpId)
        counterIndex = -1
        #print(record.reference_name,read.orgSeqName)

        if record.reference_name == read.orgSeqName:
            ## correct mapping
            counterIndex = 0  
        else:
            if orgReadSpId == mappingSpId :
                ## correct Sp but wrong mapping
                counterIndex = 1  ## map to the it org sp but map to the wrong contig
            else:
                counterIndex = 2 ## map to the wrong sp
        nhTagValue = 0
        if nhTag!= None :
            if record.is_proper_pair:
                nhTagValue = nhTag[record.qname]
            else:
                pairName =  record.qname+"_1" if record.is_read1 else   record.qname+"_2"
                nhTagValue = nhTag[pairName]
        else:
            nhTagValue = record.get_tag("NH")
        if nhTagValue == 1 :
            allCounter.unique[orgReadSpId][counterIndex]+=1
            if orgReadSpId !=mappingSpId: 
                allCounter.crossSpUnique[orgReadSpId][mappingSpId]+=1
                allCounter.crossSpTotal[orgReadSpId][mappingSpId]+=1
            continue
        #read = record.qname
        
        ## consider left and right
        pe_readId = read.readId
        if record.is_read1:
            pe_readId+= "_S1"
        if record.is_read2:
            pe_readId+= "_S2"
            
        if not pe_readId in reads  :
            reads[pe_readId] = []
            multiReadsClass[pe_readId] = []
            multiReadsOrgSp[pe_readId]= orgReadSpId
            multiReadsMappedTpSp[pe_readId] = []
            #print(record.tostring(bamFile))
            #bamRecords[read.readId] = []
            
        ## Primary Allignment
        if record.flag & 0x900 == 0 :
            allCounter.mPrimary[orgReadSpId][counterIndex]+=1
        else:
            allCounter.mSecondary[orgReadSpId][counterIndex]+=1
            allCounter.total[orgReadSpId] -=1
        reads[pe_readId].append(record)
        multiReadsClass[pe_readId].append(counterIndex)
        multiReadsMappedTpSp[pe_readId].append(mappingSpId)
        #continue
    ## Counter is  [#readsMultiMaptoOnlyOneSp,#readsMultiMaptoBothSp, #readsMultiMaptoOneWrongSp]
#    multiReadsCounter = [[0] * 3 for i in range(nSpecies)]
    for read in reads:
        classesSet = set(multiReadsClass[read])
        #print(read)
        orgSp = multiReadsOrgSp[read]
        # print(read , orgSp  , " : " , multiReadsMappedTpSp[read] , "\t " , multiReadsClass[read] )

        #print(orgSp)
        if len(classesSet ) == 1 :
            if classesSet.pop() <= 1  : #reads MultiMapto Only to One correct Sp
                allCounter.multiReads[orgSp][0]+=1
            else:
                spSet = set(multiReadsMappedTpSp[read])
                if len(spSet ) == 1 :
                    allCounter.multiReads[orgSp][2]+=1
                else:
                    allCounter.multiReads[orgSp][1]+=1
            for mappedToSp in set(multiReadsMappedTpSp[read]):
                if mappedToSp != orgSp :
                    allCounter.crossSpMulti[orgSp][mappedToSp]+=1
                    allCounter.crossSpTotal[orgSp][mappedToSp]+=1
        else:
            spSet = set(multiReadsMappedTpSp[read])
            if orgSp in spSet : ## readsMultiMaptoBothSp
                allCounter.multiReads[orgSp][1]+=1
            else:
                allCounter.multiReads[orgSp][2]+=1
            for mappedToSp in set(multiReadsMappedTpSp[read]):
                #if mappedToSp != orgSp :
                allCounter.crossSpMulti[orgSp][mappedToSp]+=1
                allCounter.crossSpTotal[orgSp][mappedToSp]+=1
    return allCounter,reads,
#%%
    


#%% Testing Paramters
## input args
#sampleFileName= "/data/bio/projects/simulation/bams/SC_to_concat1_Aligned.sortedByCoord.out.bam"
#sampleFileName= "/data/bio/projects/simulation/bams/sample_dna/cpar_calbAligned.sortedByCoord.out.bam"
#sampleFileName="/data/bio/projects/simulation/bams/sample_rna/calb_cpar_PE_trans_Aligned.sortedByCoord.out.bam"
#sampleFileName="/data/bio/projects/simulation/bams/sample_rna/calb_cpar_BWA.bam"

#From these bam file we need to calculate the following things:
#1. How many reads of CPAR map to CALB uniquely.
#2. How many read of CALB map to CALB uniquely.
#
#3. How many reads originated from CPAR are multimapped to only CPAR genome.
#4. How many reads originated from CALB are multimapped to only CALB genome.
#
#5. How many reads originated from CPAR are multimapped to both genome.
#6. How many reads originated from CALB are multimapped to both genome.
#7. Other things you think could be relevant?

#sp1InputFasta = "/data/bio/projects/simulation/bams/genome/clab.fa"
#sp2InputFasta = "/data/bio/projects/simulation/bams/genome/cpar.fa"
#
#sp1InputGTF="/data/bio/projects/simulation/bams/genome/CALB.gtf"
#sp2InputGTF="/data/bio/projects/simulation/bams/genome/CPAR.gtf"
#
#
#org1Name = "calb" # C_albicans_SC5314
#org2Name = "cpar" # C_parapsilosis_CDC317
#
#
#nSpecies = 2
#isRNA = False


def getReadCounters(args):
    logger = getLogger()
    speciesIds =   args.speciesIds # { org1Name : 0, org2Name  :1}
    
    
    logger.info("Reading Sequence Directory from Fasta files")
    
    spInputFastaFiles = args.genomes
    
    allSeqs = getSequencesPerOrganisms(spInputFastaFiles) ## testing was [sp1InputFasta,sp2InputFasta]
    seqsIndex = createSequenceIndex(allSeqs)
    seqsIndexToSeq = dict((v, k) for k,v in seqsIndex.items())
    seqToOrg = sequenceToOrganism(allSeqs)
    transcriptMap = None
    if args.simulation_type == "RNA" :
        logger.info("Reading GTF/GFF files for transcripts info ... ")
        transcriptMap = mapTranscriptToSequence( [args.annotationsGTFConcat], seqsIndex)  ## test was [sp1InputGTF, sp2InputGTF]
     
    
    counters = {}
    
    
    outputFile = open(os.path.join(args.out_dir,"report.txt"),"w")
    for rlen,layout_files in args.mappingOutputFiles.items():
        counters[rlen] = {}
        for layout,inBamFileName in layout_files.items(): 
            logger.info(f"Start counting reads for read lenghth {rlen} and ({layout}) layout:")
            
            bamFile = pysam.AlignmentFile(inBamFileName,"rb")
            ## chech if bamFile has NH tag or not , if not calc it in advance and pass it to the count method
            NHTags = checkNHTag(bamFile)
            allCounter,reads = countReads(bamFile, speciesIds , seqsIndexToSeq , seqToOrg , transcriptMap = transcriptMap , nhTag = NHTags)
            counters[rlen][layout] = allCounter
            
            outputFile.write(f"Summary Counter for lenghth {rlen} and ({layout}) layout : {inBamFileName}\n")
            allCounter.summary(outputFile)
            outputFile.write("*"*50 + "\n\n")
    outputFile.close()
    return counters


