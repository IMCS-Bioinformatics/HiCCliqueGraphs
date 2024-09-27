import csv
import json, os
from intervaltree import Interval, IntervalTree
import pandas as pd
import json, csv
import numpy as np
## Code to translate hg19->hg38 function
from pyliftover import LiftOver #https://pypi.org/project/pyliftover/
#Example usage of hg19 to hg38
# lo = LiftOver('hg19', 'hg38')
# lo.convert_coordinate('chr1', 1000000)
import re


# Directory with the current .py module
currentDir = os.path.dirname(os.path.abspath(__file__))
# Go one directory up
parentDir = os.path.dirname(currentDir)
# Navigate to the "processDifferentDB" folder
processDifferentDBDir = os.path.join(parentDir, "processDifferentDB")
# Locate the "data" directory inside "processDifferentDB"
dataDir = os.path.join(processDifferentDBDir, "data")

templ3 = { #Jung, I., Schmitt, A., et al. A compendium of promoter-centered long-range chromatin interactions in the human genome.
    "dataset": "tissuePCHiC",
    "encodeTissueMappingFn": "ENCODE_tissue-mapping.csv", #located in dataDir
    "Ufn": "encodeCcreCliqueExpressions/data/data-pvalue-0.7-fin-digraph-GRCh38.json",
    "encodeDir": f"{dataDir}/ENCODE/",
    "geneFn": f"{dataDir}/ensembl_genes_GRCh38.txt",
    "DX": 2000, #How close the genes or segments need to be for an overlap to exist
    "GTExfn": f"{dataDir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
    "GTExTissueTranslationFn": f"{dataDir}/GTEx-tissue-mapping.csv",
    #for Fantom5
    "dataFn_FANTOM5": f"{dataDir}/FANTOM5-tissue.csv",
    "translationFn": f"{dataDir}/TissuesColumns.csv", #different for 3DIV
    "resultFn": "encodeCcreCliqueExpressions/results/data-pvalue-0.7-fin-digraph-GRCh38-wEncodeBits-wExpressions-wC3.json"
}

templ2 = { #Kim, K., Jang, I., et al. 3DIV update for 2021: a comprehensive resource of 3D genome and 3D cancer genome.
    "dataset": "3DIV",
    "encodeTissueMappingFn": "ENCODE_tissue-mapping3DIV.csv", #located in dataDir
    "Ufn": "encodeCcreCliqueExpressions/data/data-pvalue-10-fin-GRCh38.json",
    "encodeDir": f"{dataDir}/ENCODE/",
    "geneFn": f"{dataDir}/ensembl_genes_GRCh38.txt",
    "DX": 2000, #How close the genes or segments need to be for an overlap to exist
    "GTExfn": f"{dataDir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
    "GTExTissueTranslationFn": f"{dataDir}/GTEx-tissue-mapping3DIV.csv",
    #for Fantom5
    "dataFn_FANTOM5": f"{dataDir}/FANTOM5-tissue.csv",
    "translationFn": f"{dataDir}/TissuesColumns3DIV.csv", #different for 3DIV
    "resultFn": "encodeCcreCliqueExpressions/results/data-pvalue-10-fin-digraph-GRCh38-wEncodeBits-wExpressions-wC3.json"
}
## !!!! Please choose the desired dataset below: !!!! ##
templ = templ3 


################################################################################################
################################################################################################
################################################################################################

geneFn = templ["geneFn"]
DX = templ["DX"]
GTExfn=templ["GTExfn"]
GTExTissueTranslationFn=templ["GTExTissueTranslationFn"]
dataFn = templ["dataFn_FANTOM5"]
translationFn=templ["translationFn"]


def addENCODE(templ):
    #Adds cCRE properties to Universal file
    #Finds cCRE properties in 
    Ufn = templ["Ufn"]
    encodeTissueMappingFn=templ["encodeTissueMappingFn"]


    #---------------------------------------------------------------------
    #Utility function that takes a file, a function that returns [ch,A,B,property] 
    #Also takes a list of segemnts.
    #Returns another list with corresponding porperties for each segment.
    #Afterwards, user may do whatever they want. Either save as a fina list or call this for each tissue type and aggregate the results
    #Another param - countabl - bool. If true, just count the rows
    #And the chromosome
    def getPropertiesForListOfSegments(fileName, segmentList, F, ch, countable=False):
        if fileName=="": 
            print ("File name empty")
            return [[] for _ in range(len(segmentList))]
        tree = IntervalTree()
        for i, segment in enumerate(segmentList):
            if len(segment)>0: tree[segment[0]:segment[1]]=i
        
        def getOverlap(min1, max1, min2, max2):
            #given 2 segments, return their overlap
            return max(0, min(max1, max2) - max(min1, min2))
        
        ##read the encode datafile
        df = pd.read_csv(fileName, compression='gzip', header=None, sep='\t', dtype=str)
        data_as_list = df.values.tolist()

        tmm = [int(r[2])-int(r[1]) for r in data_as_list]

        if countable: 
            result = [0 for _ in range(len(segmentList))] #will count number of rows in file for each segment
        else:
            result = [set() for _ in range(len(segmentList))] #list of list for each segments

        for index, row in enumerate(data_as_list):
            tmp = F(row)
            if len(tmp)==0: continue
            [chrom,A,B,props] = tmp
            if chrom!=ch: continue
            if B<A: A,B=B,A
            overlappingSegments = tree[A:B]

            if len(overlappingSegments)>0:
                overlappingSegments=list(overlappingSegments)
                for ovseg in overlappingSegments:
                    lenFeature = B-A
                    lenSegment = ovseg[1]-ovseg[0]
                    overlap = getOverlap(A,B,ovseg[0],ovseg[1])
                    if overlap<0.2*(min(lenFeature,lenSegment)):
                        continue #This overlap is too small
                    
                    if countable: result[ovseg[2]]+=1
                    else:
                        result[ovseg[2]].update(props)
        return result


    #
    def getFileName(tissue, tissueToIDMappingFn, dir):
        df = pd.read_csv(dir+tissueToIDMappingFn, usecols=["Code", "File_ID"], sep='\s*,\s*', engine='python')
        dataRow = df[df['Code'] == tissue]
        if dataRow.empty:
            print(f"No data file found with {tissue} tissue")
            return ""
        else:
            ffnn = dataRow["File_ID"].values[0]
        #If empty, then return "", processed later, all props will be empty for the tissue
        if type(ffnn)==float: return ""
        
        ffnn = dir+ffnn+".bed.gz"
        return ffnn

    def encodeGetter(row):
        #col_names = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "featureType", "classification"]
        return ([row[0], int(row[1]), int(row[2]), row[9].split(',')]) 

    #Adds information with encode cCRE props
    def addEncode(data, tfn=f"{encodeTissueMappingFn}", dir=templ["encodeDir"]):
        U=data
        possibleValues = set()

        for ch in data["chrNames"]:
            data["chrValues"][ch]["encodeStates"]=[dict() for _ in range(len(data["chrValues"][ch]["segments"]))]

            
            for tis in U["tissueIDs"]:
                print(f"Adding ENCODE for {ch}, {tis}")
                segments = U["chrValues"][ch]["segments"]
                segs = segments
                RR = getPropertiesForListOfSegments(getFileName(tis,tfn,dir), segs, F=encodeGetter, ch=ch, countable=False)
                for i,segmentSetOfProps in enumerate(RR):
                    if len(segmentSetOfProps)==0: continue
                    data["chrValues"][ch]["encodeStates"][i][tis]=list(segmentSetOfProps)
                    for el in segmentSetOfProps: possibleValues.add(el)
        data["possibleEncodes"] = list(possibleValues)
        return data


    #Read Universal
    with open(Ufn, 'r') as f:
        U = json.load(f)

    addEncode(U)

    return U

U = addENCODE(templ)
print(f"Finished adding ENCODE cCRE properties")

#######################################################################


##Adding genes
##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def addGenes(U, geneFn = geneFn):
    #First some utility functions, then the loop over all chrs to add genes 
    #---------------------------------------------------------------------
    #Utility function that takes a file with name fileName, a function F that returns [ch,A,B,property] 
    #Also takes a list of segemnts segmentList.
    #Returns another list with corresponding properties for each segment.
    #Afterwards, user may do whatever they want. Either save as a final list or call this for each tissue type and aggregate the results
    #Another param - countable - bool. If true, just count the number of rows that overlap with each of the segments
    #And the chromosome ch
    #Adjusted for ensembl_genes_GRCh38.txt file processing, works only with this file
    def getPropertiesForListOfSegments(fileName, segmentList, F, ch, countable=False, dx=DX):
        if fileName=="": 
            return [[] for _ in range(len(segmentList))]
        
        tree = IntervalTree() #The graph segments are stored in this tree to find other intervals that overlap with segments fast
        for i, segment in enumerate(segmentList):
            if len(segment)>0: tree[segment[0]-dx:segment[1]+dx]=i
        
        def getOverlap(min1, max1, min2, max2):
            #given 2 segments, return their overlap
            return max(0, min(max1, max2) - max(min1, min2))
        
        def getBroadGeneType(type):
            #returns an object, e.g. {'protein coding': True, 'RNA': false, 'pseudogene':False}
            return {'protein coding': ('protein_coding' in type), 'RNA': ('RNA' in type), 'pseudogene':('pseudo' in type)}
        
        ##read the ensembl genes datafile
        #Its filename is in fileName argument
        df = pd.read_csv(fileName, header=1, sep=',', dtype=str)
        data_as_list = df.values.tolist()

        tmm = [int(r[4])-int(r[3]) for r in data_as_list if r[2].isdigit()]
        iii=1
        if countable: 
            result = [0 for _ in range(len(segmentList))] #will count number of rows in file for each segment
        else:
            result = [set() for _ in range(len(segmentList))] #list of sets for each segment
        geneIDToGeneInfo = {} #All possible information about a gene is mapped to its ID. Only those genes are recorded that overlap with at least one segment

        for index, row in enumerate(data_as_list):
            tmp = F(row) #read the next row from properties file (e.g. file with all existing ensembl genes) and extract the useful information - ch, start, end, otherData
            if len(tmp)==0: continue
            [chrom,A,B,geneID] = tmp
            if (chrom!=ch) and (str(chrom)!=str(ch)) and (f"chr{chrom}"!=ch): continue # if chromosomes are different from what is currently processed, go to next gene
            if B<A: A,B=B,A
            overlappingSegments = tree[A:B] #all segments that overlap with the current gene

            if len(overlappingSegments)>0:
                iii=9
                overlappingSegments=list(overlappingSegments)
                for ovseg in overlappingSegments:
                    lenFeature = B-A #length of a feature from propertiy file. E.g., length of a gene
                    lenSegment = ovseg[1]-ovseg[0] #Hi-C segment length
                    overlap = getOverlap(A,B,ovseg[0],ovseg[1])
                    # if overlap<0.2*(min(lenFeature,lenSegment)):
                    #     continue #This overlap is too small #removed 
                    
                    if countable: result[ovseg[2]]+=1
                    else:
                        result[ovseg[2]].add(geneID) #map the gene to the segment that it overlaps
                        #..and store the full information about the gene to geneIDToGeneInfo
                        if geneID in geneIDToGeneInfo:
                            if geneIDToGeneInfo[geneID]["length"]!=abs(int(row[4])-int(row[3])):
                                raise Exception("Two genes with the same ID have different lengths, sanity check failed")
                            continue #all good, the record is already saved before
                        geneIDToGeneInfo[geneID] = {
                            #row==[ "Gene stable ID","Gene stable ID version","Chromosome/scaffold name",
                                    # "Gene start (bp)","Gene end (bp)","Strand","Gene name","Gene type"]
                            "ID": row[0],
                            "name": str(row[6]),
                            "start": int(row[3]),
                            "end": int(row[4]),
                            "length": abs(int(row[4])-int(row[3])),
                            "strand": row[5],
                            "type": row[7],
                            "broadType": getBroadGeneType(row[7]),
                        }
                        if geneIDToGeneInfo[geneID]["name"]=="nan": #some genes have no names
                            geneIDToGeneInfo[geneID]["name"]=""
        return [result, geneIDToGeneInfo]

    def ensemblGeneGetter(row):
        #For parsing a line from ensembl_genes_GRCh38.txt
        #col_names = ["Gene stable ID","Gene stable ID version","Chromosome/scaffold name",
        # "Gene start (bp)","Gene end (bp)","Strand","Gene name","Gene type"]
        return ([str(row[2]), int(row[3]), int(row[4]), row[0]]) 

    ###Start the main loop to add genes
    for ch in U["chrNames"]:
        print("Adding ensembl genes for ", ch)
        segments = U["chrValues"][ch]["segments"]
        U["chrValues"][ch]["segmentGenes"]=[[] for _ in range(len(segments))] #to each segment genes are mapped, using indeces from chromosomeGenes list

        iii=0
        [res, geneIDToGeneInfo] = getPropertiesForListOfSegments(geneFn, segments, ensemblGeneGetter, ch)
        iii+=1
        listOfChromosomeGenes = sorted([gene for gene in geneIDToGeneInfo.values()], key=lambda x:(x["start"], x["end"]))
        U["chrValues"][ch]["chromosomeGenes"]=listOfChromosomeGenes #each list element is a dict with data about this gene
        iii+=1
        geneIDToIndex = {gene["ID"]:i for i,gene in enumerate(listOfChromosomeGenes)}
        
        for i,setOfSegmentGenes in enumerate(res):
            for gene in setOfSegmentGenes: #for each gene that overlaps the current segment
                U["chrValues"][ch]["segmentGenes"][i].append(geneIDToIndex[gene]) #Append its index in the all genes array for the chromosome to the list
    return U
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#Add GTEx to genes
#U["chrValues"][ch]["chromosomeGenes"] is a list of genes in the chromosome.
#Each gene record is a dict. The goal now is to add a new entry to each gene record with data on its GTEx expression

def addGeneGTExData(U, GTExfn=GTExfn, GTExTissueTranslationFn=GTExTissueTranslationFn):

    #read GTEx file to G
    with open(GTExfn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        G = list(spamreader)
        G = G[2:] #remove first 2 lines
        for i,row in enumerate(G):
            G[i][0] = row[0].split('.')[0]
        headerG = G.pop(0) #remove the header

        ensIDToRowInG = {row[0]:i for i,row in enumerate(G)}

    
    #read file with tissue to column name
    with open(GTExTissueTranslationFn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        sprl = list(spamreader)
        tisToTitleRows = {row[0].strip():row[1] for row in sprl}
        titleRowsToTis = {row[1]:row[0].strip() for row in sprl}
        tisToHeaderIndex = {tis: headerG.index(fullTissueName) for tis,fullTissueName in tisToTitleRows.items() if fullTissueName in headerG}

    for ch in U["chrNames"]:
        for i,gene in enumerate(U["chrValues"][ch]["chromosomeGenes"]):
            geneID = gene['ID']
            U["chrValues"][ch]["chromosomeGenes"][i]["GTEx"]={}
            if geneID not in ensIDToRowInG:
                continue
            gtexRecord = G[ensIDToRowInG[geneID]]
            for tis in U['tissueIDs']:
                if tis not in tisToHeaderIndex:
                    continue
                U["chrValues"][ch]["chromosomeGenes"][i]["GTEx"][tis]=float(gtexRecord[tisToHeaderIndex[tis]])
                iii=0
        iii=0 #end looping over one chromosome's genes
    #end loop over chromosomes
    return U


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#Add fantom5 profiles
def addFantoms(U, dataFn=dataFn, translationFn=translationFn, dx=DX):

    #Read data file with data on Fantom5
    with open(dataFn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        L = list(spamreader)

    #Read translation file tissue name to column name
    with open(translationFn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        T = list(spamreader)

    iii=0


    
    def parse_genomic_string(s):
        """
        Parses a genomic string and returns its components.

        Args:
        s (str): A genomic string in the format 'genome::chromosome:start..end,strand;extra_info'

        Returns:
        dict: A dictionary containing the parsed components.
        """
        pattern = r'(\w+)::(chr[\w\d]+):(\d+)\.\.(\d+),([+-]);(.+)'
        match = re.match(pattern, s)

        if not match:
            raise ValueError("String does not match the expected format")

        # Example usage
        # genomic_string = "hg19::chr1:565266..565278,+;hg_3.1"
        # parsed_data = parse_genomic_string(genomic_string)
        # print(parsed_data)
        return {
            'genome': match.group(1),
            'chromosome': match.group(2),
            'start': int(match.group(3)),
            'end': int(match.group(4)),
            'strand': match.group(5),
            'extra_info': match.group(6)
        }
    #process all chromosomes

    for ch in U["chrNames"]:
        print("Adding FANTOM5 for ", ch)
        #Finding bitmaps for segments
        allSegments = U["chrValues"][ch]["segments"]
        allLinks = U["chrValues"][ch]["links"]
        
        curChrDataRows=[]
        for i,dataRow in enumerate(L):
            if i==0: continue
            info=parse_genomic_string(dataRow[0])
            if info["chromosome"]!=ch and False:
                continue
            else:
                info["realIndex"]=i #row nr in data file
                info["expressionInTissues"]={}
                for tissue in U["tissueBits"].keys():
                    columnIDs = [j for j in range(len(T)) if T[j][0].strip()==tissue]
                    if len(columnIDs)==0: continue
                    fantomExpressionInCurrentTissue = max([float(dataRow[colID]) for colID in columnIDs])
                    info["expressionInTissues"][tissue] = fantomExpressionInCurrentTissue
                info["description"] = dataRow[1]
                info["uniprot id"] = dataRow[5]
                # info["association with transcript"] = dataRow[2]
                curChrDataRows.append(info)
                
        curChrDataRows = sorted(curChrDataRows, key=lambda x: (x["start"], x["end"]))
        iii=0
        #Translate the start and end coordinates to GRCh38
        ## Code to translate hg19->hg38 function
        # from pyliftover import LiftOver #https://pypi.org/project/pyliftover/
        # #Example usage of hg19 to hg38
        # lo = LiftOver('hg19', 'hg38')
        # lo.convert_coordinate('chr1', 1000000)
        lo = LiftOver('hg19', 'hg38')
        grch38Fantoms = []
        for j,record in enumerate(curChrDataRows):
            start,end = lo.convert_coordinate(ch, record['start']), lo.convert_coordinate(ch, record['end'])
            #If some of the endpoints have no corresponding locus in hg38, throw it away
            if len(start)==0 or len(end)==0:
                # print("failed to translate a locus")
                continue
            
            #If segments change the chromosome, throw it away
            if start[0][0]!=ch or end[0][0]!=ch:
                print("bad chr")
                continue
                
            
            #Nomal case, assign new loci
            [start, end]=sorted([start[0][1],end[0][1]])
            record["start"], record["end"] = [start, end]
            record['genome']="GRCh38"
            grch38Fantoms.append(record)
        iii=0

        tree = IntervalTree() #The graph segments are stored in this tree to find other intervals that overlap with segments fast
        for i, segment in enumerate(allSegments):
            if len(segment)>0: tree[segment[0]-dx:segment[1]+dx]=i

        fantomsOfSegments = [set() for _ in range(len(allSegments))] #Each set corresponds to those fantoms that are overlapping it. Each set contains record indeces from grch38Fantoms list
        # chromosomeFantoms = {}#realIndex: record

        for index, fantomRecord in enumerate(grch38Fantoms):
            chrom= fantomRecord["chromosome"]
            if (chrom!=ch) and (str(chrom)!=str(ch)) and (f"chr{chrom}"!=ch): continue # if chromosomes are different from what is currently processed, go to next gene
            A,B = fantomRecord["start"], fantomRecord["end"]
            overlappingSegments = tree[A:B] #all segments that overlap with the current gene

            if len(overlappingSegments)>0:
                iii=9
                overlappingSegments=list(overlappingSegments)
                for ovseg in overlappingSegments:
                    iii=0
                    segmentID = ovseg[2]
                    fantomsOfSegments[segmentID].add(index)

        U["chrValues"][ch]["segmentFantoms"] = [[] for _ in range(len(allSegments))]
        U["chrValues"][ch]["chromosomeFantoms"] = []

        fantomsInAtLeastOneSegment  = sorted(list(set.union(*fantomsOfSegments)))
        indexInAllFantomsToIndexInActiveFantoms = {}
        activeFantoms = []
        for i, fantomIndex in enumerate(fantomsInAtLeastOneSegment):
            fantomRecord = grch38Fantoms[fantomIndex]
            activeFantoms.append(fantomRecord)
            indexInAllFantomsToIndexInActiveFantoms[fantomIndex]=i

        U["chrValues"][ch]["chromosomeFantoms"] = activeFantoms
        keysToRemove = ["genome", "realIndex", "chromosome", "realIndex", "uniprot id"]
        U["chrValues"][ch]["chromosomeFantoms"] = [{k: v for k, v in d.items() if k not in keysToRemove} for d in U["chrValues"][ch]["chromosomeFantoms"]]

        for seg,fantomSet in enumerate(fantomsOfSegments):
            fantomList = []
            for el in fantomSet:
                fantomList.append(indexInAllFantomsToIndexInActiveFantoms[el])
            U["chrValues"][ch]["segmentFantoms"][seg]=fantomList

    return U

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

U = addGenes(U)
print("added genes")
U = addGeneGTExData(U)
print("added GTEx")
U = addFantoms(U)
print("Added FANTOM5")

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#This part will combine ENCODE annotations into 4 groups and give them bitmaps instead of lists of str
#Groups are:
# (PLS or pELS -->promoter), 
# (dELS --> distal enhancer), 
# (CTCF-bound or CTCF-only --> CTCF), 
# (DNase-only)
def addEncodeBits(U):

    encodeToBitmap = {"PLS":1,
                        "pELS": 1,
                        "dELS":2,
                        "CTCF-bound":4,
                        "CTCF-only":4,
                        "DNase-only":8,
                        }
    encodeBits = {
        "promoter":1,
        "distal enhancer":2,
        "CTCF":4,
        "DNase-only":8,
    }

    def getBitmapForListOfEncodes(listOfEncodes):
        curBit=0
        for enc in listOfEncodes:
            curBit|=encodeToBitmap.get(enc, 0)
        return curBit

    tissueWithEncode = set()
    for ch in U["chrNames"]:
        for i,encDict in enumerate(U["chrValues"][ch]["encodeStates"]):
            if len(encDict)==0: continue
            newDict = {}
            for tis,listOfEncodes in encDict.items():
                newBit = getBitmapForListOfEncodes(listOfEncodes)
                if newBit>0:
                    newDict[tis]=newBit
                    tissueWithEncode.add(tis)
            U["chrValues"][ch]["encodeStates"][i]=newDict
            iii=0
    U["encodeBits"] = encodeBits
    del U["possibleEncodes"]

    #Find tissues that have info on fantom expressions and on genes and on GTEx
    tissueWithFantoms = set()
    for ch in U["chrNames"]:
        for i,fanDict in enumerate(U["chrValues"][ch]["chromosomeFantoms"]):
            curTissues = fanDict["expressionInTissues"].keys()
            curTissues = set(curTissues)
            tissueWithFantoms.update(curTissues)

    tissuesWithGTEx = set()
    for ch in U["chrNames"]:
        for i,geneDict in enumerate(U["chrValues"][ch]["chromosomeGenes"]):
            curTissues = geneDict["GTEx"].keys()
            curTissues = set(curTissues)
            tissuesWithGTEx.update(curTissues)



    U["dataAviabilityForTissues"]={
        "encodeStates": sorted(list(tissueWithEncode)),
        "GTEx": sorted(list(tissuesWithGTEx)),
        "fantoms": sorted(list(tissueWithFantoms)),
    }

    return U



##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
U = addEncodeBits(U)
print(f"Finished adding bits for encode properties")

### Add and calculate cliques of size 3
def addC3(U):
    def makeAdjAndBitsOfLink(links):
        # creates an adj structure for all links and creates a dictionary bitsOfLink
        # bitsOfLink[(A,B)] == the bitmap of the link
        adj = dict()
        bitsOfLink = dict()
        for link in links:
            [A, B, bitmap] = link
            if (B<=A): 
                tmm=A
                A=B
                B=tmm
            if A not in adj: adj[A] = set()
            adj[A].add(B)
            if B not in adj: adj[B] = set()
            adj[B].add(A)
            bitsOfLink[(A,B)] = bitmap
        return [adj, bitsOfLink]
        
    def getC3(links):
        # finds all triangles in current chr. At least one tissue must be shared among each link of the triangle
        # result is a list of tuples (A, B, C, bit), where A, B, C are sorted asc. and where bit is the max bitmap shared among all links (AB)&(AC)&(BC)
        # only those triangles are kept that have at least one common tissue in all 3 links
        # uses links that may be directed
        # links argumentthat is used here has list of links in form [A, B, bitmap] for one chromosome
        
        [adj, bitsOfLinks] = makeAdjAndBitsOfLink(links) 
        #adj[A] == set of all other vertices that are adjacent to A
        #bitsOfLink[(A,B)] == the bitmap of the link (A,B)
        cliques = set() #{(A,B,C,bit), (), ()}
        
        # Algorithm to find C3: take vertice A. setOfBs has all its direct neighbors (B>A to get A<B<C in result and get rid of duplicates)
        # Each node in setOfBs can be on a C3. If A-B is in some C3, then there exists a C that is both adjacent to B and to A.
        # Nodes that are adjacent to A are stored in setOfBs. Nodes that are adjacent to current B are stored in setOfCs.
        # Taking intersection - if it is non-empty, for each C, A-B-C is a triangle (clique)
        for A in adj.keys():
            setOfBs = adj[A]
            setOfBs = set(el for el in setOfBs if el>A)
            for B in setOfBs:
                setOfCs = adj[B]
                setOfCs = set(el for el in setOfCs if el>B)
                
                goodCs = setOfCs.intersection(setOfBs)
                for C in goodCs:
                    [A,B,C] = sorted([A,B,C])
                    bits = ( bitsOfLinks[(A,B)] & bitsOfLinks[(B,C)] & bitsOfLinks[(A,C)])
                    if bits>0: #why 0? if want all CL3 not looking at tissue shared, can ignore this if
                        cliques.add((A,B,C,bits))
        del adj
        del bitsOfLinks
        return cliques  

    for ch in U["chrNames"]:
        print("Calculating C3 for ", ch)
        links = [link[:3] for link in U["chrValues"][ch]["links"]]
        iii=0
        cliques = getC3(links)
        iii+=1
        listOfCliques = sorted(list(cliques), key=lambda x: (x[0], x[1], x[2]))
        U["chrValues"][ch]["cliques"] = listOfCliques
    return U


##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
U = addC3(U)

with open(templ["resultFn"], 'w') as file:
    json.dump(U, file)
