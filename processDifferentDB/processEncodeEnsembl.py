#This module adds encode and ensembl data to the Hi-C graph nodes
#Data files that are processed here will be taken from:
    # the /graphFiles directory to get the Hi-C graphs 
    # the /data directory to get the encode and ensembl data for different tissue types
    # The tissue mapping to file names is also stroed in /data

import os
# Get the directory of the current script
current_dir = os.path.dirname(os.path.abspath(__file__))



# Directory with DB data
target_dir_name = "data"
dataDir = os.path.join(current_dir, target_dir_name)

templ2 = { #Kim, K., Jang, I., et al. 3DIV update for 2021: a comprehensive resource of 3D genome and 3D cancer genome.
    "dataset": "3DIV",
    "dir_name_with_graphs": "graphFiles3DIV",
    "dataFn": f"{dataDir}/FANTOM5-tissue.csv",
    "translationFn": f"{dataDir}/TissuesColumns3DIV.csv",
    "encodeTissueMappingFn": "ENCODE_tissue-mapping3DIV.csv", #located in dataDir
    "GTExTissueMappingFn":"GTEx-tissue-mapping3DIV.csv", #in dataDir
    "ENCODE_RNA_tissue_mappingFn": "ENCODE_RNA_tissue-mapping3DIV.csv", #in dataDir/ENCODE
    "CALC_CENTR": False,
    "CHRS": ["chrX"]+[f"chr{i}" for i in range(1,23,1)]
}
templ3 = { #Jung, I., Schmitt, A., et al. A compendium of promoter-centered long-range chromatin interactions in the human genome.
    "dataset": "tissuePCHiC",
    "dir_name_with_graphs": "graphFilesTissuePCHiC",
    "dataFn": f"{dataDir}/FANTOM5-tissue.csv",
    "translationFn": f"{dataDir}/TissuesColumns.csv",
    "encodeTissueMappingFn": "ENCODE_tissue-mapping.csv", #located in dataDir
    "GTExTissueMappingFn":"GTEx-tissue-mapping.csv", #in dataDir
    "ENCODE_RNA_tissue_mappingFn": "ENCODE_RNA_tissue-mapping.csv", #in dataDir/ENCODE
    "CALC_CENTR": False,
    "CHRS": ["chrX"]+[f"chr{i}" for i in range(1,23,1)]
}
templ = templ2 #choose one of the datasets

DS = templ["dataset"]

# Directory with Hi-C graphs
target_dir_name = templ["dir_name_with_graphs"]
udir = os.path.join(current_dir, target_dir_name)

#Adding FANTOM5 data
dataFn=templ["dataFn"]
translationFn = templ["translationFn"]

encodeTissueMappingFn=templ["encodeTissueMappingFn"]
GTExTissueMappingFn=templ["GTExTissueMappingFn"]
ENCODE_RNA_tissue_mappingFn = templ["ENCODE_RNA_tissue_mappingFn"] 

CHRS = templ["CHRS"]
CALCULATE_CENTRALITIES = templ["CALC_CENTR"] #it is slow and not necessary - it is an option to calculate centralities of nodes
#CHRS = 

import csv, json

#Read data file with data on Fantom5
with open(dataFn, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
    L = list(spamreader)

#Read translation file tissue name to column name
with open(translationFn, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    T = list(spamreader)


import re
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

from statistics import median
import numpy as np
import statistics

def get_basic_statistics(data):
    if not data:
        return "List is empty"

    # Calculate mean
    mean = statistics.mean(data)

    # Calculate quartiles
    quartiles = np.percentile(data, [25, 50, 75])

    # Calculate max
    max_value = max(data)

    return {
        "mean": mean,
        "25%": quartiles[0],
        "50%": quartiles[1],  # This is the median
        "75%": quartiles[2],
        "max": max_value
    }

fullStats = []

import matplotlib.pyplot as plt

#These functions process the data
def getFantomkExpressionValue(rowId, tissue):
    #returns fantom expression value from rowId of L, and tissue
    #T is a translation file, read before
    if rowId<0:
        raise Exception(f"rowId can not be negative {rowId}")
    columnIDs = [j for j in range(len(T)) if T[j][0].strip()==tissue]
    values = []
    for colID in columnIDs:
        values.append(float(L[rowId][colID]))
    return max(values)

def getFantomName(rowId):
    if rowId<0:
        raise Exception(f"rowId can not be negative {rowId}")
    columnID = 1 #take 2nd column
    return L[rowId][columnID] #e.g. "p1@ACTB"

def getMinDistancesForAllSegments(segments, fantoms):
    i=0
    j=0 # currentHypothesisForClosestFantom
    def skelas(a,b): #overlap?
        return a[0]<=b[1] and b[0]<=a[1]  #!!!
    def isOnTheLeft(star, makonis): #star is a fantom segment, makonis is a segment. 
        return star[1]<makonis[0]
    fantoms.append([10e14, 10e15+10,-1]) #append a dummy fantom that is definetly not the closest one to any segment
    res = [10e10 for _ in range(len(segments))] #for each segment and for each tissue
    fantomIndeces = [[] for _ in range(len(segments))] #for each segment, find the name of the fantom, use 2nd column of the data file
    while (i<len(segments)) :
        makonis = segments[i]
        star = fantoms[j]

        if skelas(makonis, star):
            
            if res[i]>0: #ja tur jau bija nenulle
                fantomIndeces[i]=[] #ja pirms tam bija kaut kas aizpildīts un tgd ieraugām kaut ko, kas šķeļas, to apnuļļo
            res[i]=0
            fantomIndeces[i].append(star[2])

            J=j
            STAR=star
            while(J<len(fantoms)-1):
                J+=1
                STAR=fantoms[J]
                if skelas(makonis, STAR):
                    fantomIndeces[i].append(STAR[2])
                else: 
                    break
            i+=1
            continue
        
        if (isOnTheLeft(star, makonis)):
            res[i]=min(res[i],makonis[0]-star[1])
            fantomIndeces[i]=[star[2]]
            j+=1
            continue
        dist = star[0]-makonis[1]
        
        if j>0 and makonis[0]-fantoms[j-1][1] < dist :
            #previous fantom was closer
            dist=makonis[0]-fantoms[j-1][1]
            fantomIndeces[i]=[fantoms[j-1][2]]
        else:
            fantomIndeces[i]=[fantoms[j][2]]
        res[i] = dist

        i+=1
    fantoms.pop(-1) #remove the dummy element from the back. ~10:10
    return [res,fantomIndeces]



#process all chromosomes
for ch in CHRS:
    fn = f"{udir}/{ch}.json"
    with open(fn, 'r') as f:
        U = json.load(f)
    #Finding bitmaps for segments
    allSegments = U["chrValues"][ch]["segments"]
    allLinks = U["chrValues"][ch]["links"]
    
    segmentBitmaps = [0 for _ in range(len(allSegments))]
    for [a,b,bit] in allLinks:
        segmentBitmaps[a]|=bit
        segmentBitmaps[b]|=bit
    
    
    curChrDataRows=[]
    for i,dataRow in enumerate(L):
        if i==0: continue
        info=parse_genomic_string(dataRow[0])
        if info["chromosome"]!=ch:
            continue
        else:
            info["realIndex"]=i #row nr in data file
            info["tissueBitmap"]=0
            resFantomTissueBitmap = 0
            for tissue in U["tissueData"]["tissueBits"].keys():
                tisBit = U["tissueData"]["tissueBits"][tissue]
                columnIDs = [j for j in range(len(T)) if T[j][0].strip()==tissue]
                for colID in columnIDs:
                    if float(dataRow[colID])>0 or True: #Changed here - now expression=0 is a valid result as well
                        info["tissueBitmap"]|=tisBit
            curChrDataRows.append(info)
            
    curChrDataRows = sorted(curChrDataRows, key=lambda x: (x["start"], x["end"]))

    fantomResults = [{} for _ in range(len(allSegments))]
    fantomValues = [{} for _ in range(len(allSegments))]
    fantomNames = [{} for _ in range(len(allSegments))]

    for tissue in U["tissueData"]["tissueBits"].keys():
        tisBit = U["tissueData"]["tissueBits"][tissue]
        tissueSegments = [allSegments[i] for i in range(len(allSegments)) if ((segmentBitmaps[i]&tisBit)==tisBit)]
        tissueSegmentInds = [i for i in range(len(allSegments)) if ((segmentBitmaps[i]&tisBit)==tisBit)]
    
        fantomSegments = \
            [[curChrDataRows[i]["start"],curChrDataRows[i]["end"],curChrDataRows[i]["realIndex"]] for i in range(len(curChrDataRows)) if ((curChrDataRows[i]["tissueBitmap"]&tisBit)==tisBit)]
        
        #fantomSegmentInds = [i for i in range(len(curChrDataRows)) if ((curChrDataRows[i]["tissueBitmap"]&tisBit)==tisBit)]

        [distToF, fantomTissueIndeces] = getMinDistancesForAllSegments(tissueSegments, fantomSegments)
        iii=0

        
        for k,dist in enumerate(distToF):
            #k is index, dist is the distance from segment to the closest fantom
            #store results to fantomResults list of dicts
            segIndInAllSegments = tissueSegmentInds[k] #distToF only has data on segments that are in the current tissue. Tranlating the index to index in allSegments 
            fantomResults[segIndInAllSegments][tissue] = distToF[k]
            fantomValues[segIndInAllSegments][tissue] = max([getFantomkExpressionValue(ind, tissue) for ind in fantomTissueIndeces[k]])  #bestExprs[k] 
            fantomNames[segIndInAllSegments][tissue] = [getFantomName(ind) for ind in fantomTissueIndeces[k]]#
            iii=9
            # if k>=len(curChrDataRows):
            #     iiierr=1
            # actualFantom = curChrDataRows[k] #might be used later
            # fantomRow = L[actualFantom["realIndex"]] # the records for the current fantom
            # iii=0
        iii=0 #end processing one tissue

    iii=0 #end processing all tissues
    U["chrValues"][ch]["closestFantoms"] = fantomResults
    U["chrValues"][ch]["segmentFantoms"] = fantomValues
    U["chrValues"][ch]["fantomNames"] = fantomNames

    U["VersionInfo"] = {"Generated using": "processEncodeEnsembl.py",
                              "inputFile": "",
                              "inputFileMadeBy": "",
                              "purpose": "File with hg19 properties (FANTOM5 data) added",
                              "version": "05-02-24-1",
                              }
    #######################################################################################
    ffnn = f"{udir}/rez/wfantoms05-{ch}.json"
    with open(ffnn, 'w', encoding ='utf8') as json_file: 
        json.dump(U, json_file, ensure_ascii = False) 
    print(ffnn)

#at this point fantom data is added to all chr files and new temporary files are created
#now the files will be translated to hg38 and ensembl properties (in hg38) will be added


from copy import deepcopy
from intervaltree import Interval, IntervalTree
import pandas as pd
import math

## Code to translate hg19->hg38 function
from pyliftover import LiftOver #https://pypi.org/project/pyliftover/
#Example usage of hg19 to hg38
lo = LiftOver('hg19', 'hg38')
lo.convert_coordinate('chr1', 1000000)

def translateGeneFile(data):
    #Input- contents of a file, json object, from readGeneFile function
    #Output - another json object with hg38 coordinates for segments
    ch = data["chrNames"][0]
    data["tissueData"]["tissueIDs"] = sorted(list(data["tissueData"]["tissueBits"].keys()))
    newData = deepcopy(data)

    noCorresponding = 0 #For statistic of hg19 to fg38 conversion
    otherChr = 0
    savableChr=0
    strand=0
    segments = newData["chrValues"][ch]["segments"]
    for i, seg in enumerate(segments):
        newA = lo.convert_coordinate(ch, seg[0])
        newB = lo.convert_coordinate(ch, seg[1])
 
        #If some of the endpoints have no corresponding locus in hg38, throw it away
        if len(newA)==0 or len(newB)==0:
            noCorresponding+=1
            segments[i]=[]
            continue
        
        #If segments change the chromosome, throw it away
        if newA[0][0]!=ch or newB[0][0]!=ch:
            otherChr+=1
            if newA[0][0] in data["chrNames"] and newB[0][0] in data["chrNames"] and newA[0][0]==newB[0][0]:
                #This means that both endpoints went to another neat chromsome. Could be recovered and moved to another chr ifd there are too many of theese
                savableChr+=1
            segments[i]=[]
            continue
            
        #if either segment is now on some other strand, i.e. not + strand
        if newA[0][2]!='+' or newB[0][2]!='+':
            strand+=1
            segments[i]=[]
            continue
        
        #Nomal case, assign new loci
        segments[i]=sorted([newA[0][1],newB[0][1]])

        #Check if those are int
        if type(segments[i][0])!=type(1) or type(segments[i][1])!=type(1):
            raise Exception(f"Something wrong with translation hg19 to hg38, {segments[i]}")
    print(f"hg19 to hg38 done: {noCorresponding=}, {otherChr=}, {savableChr=}, {strand=} from {len(segments)=}", ch)

    #Some segments after translation to hg38 became invalid ([]).
    #Therefore some links became invalid too.
    #I will simply remove all invalid links - thats the easy option
    #Alternative - remove [] segments and update link segment indeces - more difficult because there are aribitrarely many lists that might use segment indeces, e.g. genes or states
    segments = newData["chrValues"][ch]["segments"]
    badSegmentInds = [i for i in range(len(segments)) if segments[i]==[]]
    badSegmentInds = set(badSegmentInds)
    newData["chrValues"][ch]["links"] = [link for link in newData["chrValues"][ch]["links"] if link[0] not in badSegmentInds and link[1] not in badSegmentInds]
    
    return newData

##validifyUniversal(data) - update segemnt indeces and the list of segments
def validifyUniversal(data):
    O = data["chrValues"][data["chrNames"][0]] #pointer
    #Segments now have new indeces
    validSegmentInds = [i for i,seg in enumerate(O["segments"]) if len(seg)>0]
    validSegmentSet = set(validSegmentInds)

    local_oldIndexToNewIndex = {i: validSegmentInds.index(i) if len(seg)>0 else -1 for i,seg in enumerate(O["segments"])  }
    def oldToNewIndex(oldIndex):
        return local_oldIndexToNewIndex[oldIndex]
    
    newSegmentList = [seg for i,seg in enumerate(O["segments"]) if len(seg)>0]
    O["segments"] = newSegmentList

    #Updating links
    newLinks = []
    for link in O["links"]:
        oldA,oldB = link[0],link[1]
        if oldA in validSegmentSet and oldB in validSegmentSet:
            newLink = [el for el in link] #copy of th elink
            newLink[0] = oldToNewIndex(oldA)
            newLink[1] = oldToNewIndex(oldB)
            newLinks.append(newLink)
    O["links"] = newLinks
    #Updating segmentGenes
    O["segmentGenes"] = [el for i,el in enumerate(O["segmentGenes"]) if i in validSegmentSet]
    #Updating segmentStates
    O["segmentStates"] = [el for i,el in enumerate(O["segmentStates"]) if i in validSegmentSet]

    #Updating closestFantoms
    O["closestFantoms"] = [el for i,el in enumerate(O["closestFantoms"]) if i in validSegmentSet]
    #Updating segmentFantoms
    O["segmentFantoms"] = [el for i,el in enumerate(O["segmentFantoms"]) if i in validSegmentSet]
    #Updating fantomNames
    O["fantomNames"] = [el for i,el in enumerate(O["fantomNames"]) if i in validSegmentSet]

    return data

#------------------------------------------------------------------------------------
##Process ensembl files, add everything from them
def processAllEnsembl(data, edir=f"{dataDir}/ensembl/"):
    #the edir should contain Human_Regulatory_Features_(GRCh38.p14) among other files
    ensFns = []
    for filename in os.listdir(edir):
        filepath = os.path.join(edir, filename)
        if os.path.isfile(filepath):
            ensFns.append(filename)

    ensFnToProperty = {
        'Human_miRNA_Target_Regions_(GRCh38.p14).txt': 'Feature type class', 
        'Human_Other_Regulatory_Regions_(GRCh38.p14).txt': 'Feature type class', 
        'Human_Regulatory_Features_(GRCh38.p14).txt': 'Feature type', 
        'Human_Somatic_Short_Variants_(SNPs_and_indels_excluding_flagged_variants)_(GRCh38.p14).txt': 'Variant name',
        'Human_Somatic_Structural_Variants_(GRCh38.p14).txt': 'Structural variant name',
        'Human_Structural_Variants_(GRCh38.p14).txt': 'Structural variant name',
    }
    ensToName = {
        prop: prop[6:].split('_(')[0] for prop in ensFns
    }

    #I have ensFns with filenames in edir (EnsemblDIRectory)
    #Some files have bool type, others are to be counted
    boolableFNs=['Human_miRNA_Target_Regions_(GRCh38.p14).txt', 'Human_Other_Regulatory_Regions_(GRCh38.p14).txt', 'Human_Regulatory_Features_(GRCh38.p14).txt']
    countableFNs=['Human_Somatic_Short_Variants_(SNPs_and_indels_excluding_flagged_variants)_(GRCh38.p14).txt', 'Human_Somatic_Structural_Variants_(GRCh38.p14).txt', 'Human_Structural_Variants_(GRCh38.p14).txt']

    def getOverlap(min1, max1, min2, max2):
        #given 2 segments, return their overlap
        return max(0, min(max1, max2) - max(min1, min2))

    def processFeatureFile(ensemblFn, data, tree, property, featureName="", possibleValues=set(), chrom="chr1", boolable=True):
        """
        This function processes an ensembl file row by row. 
        If overlap with any segment in our graph is significant (>=20% of the smallest), a property will be saved for that segment
        this function adds some fields to the data:
            adds to list data["boolable properties"] or to data["countable properties"]
            creates a list of list or list of int in data["chrValues"][chrom][property]
        """
        #Call e.g.     
        # #processFeatureFile(ensemblFn="C:/Users/asizo/Documents/myPrograms/dataFiles/ensembl/Human_miRNA_Target_Regions_(GRCh38.p14).txt", 
        #   data=data, 
        #   tree=tree, 
        #   property="Human_miRNA_Target_Regions", (e.g.) 
        #   featureName="Feature type class", (e.g.)
        #   possibleValues=possibleValues, chrom=chrom)
        #
        # where tree = IntervalTree()
        # for i, segment in enumerate(data["chrValues"][chrom]["segments"]):
        #     if len(segment)>0: tree[segment[0]:segment[1]]=i
        #if boolable, then each property is either 0 or 1. Else - count of occurances
        # possibleValues is a set()
        # feature name is exactly as in the ensembl file - column name to be processed
        ##########################################################################
        if boolable:
            #Each property is either 0 or 1 for each segment
            #Each segment receives a list of properties
            if "boolable properties" not in data:
                data["boolable properties"]=[]
            data["boolable properties"].append(property)
        else:
            #Number of existing properties will be counted. 
            #No matter what's the property name, we are interested in the count only
            if "countable properties" not in data:
                data["countable properties"]=[]
            data["countable properties"].append(property) 
        
        if boolable:
            data["chrValues"][chrom][property] = [[] for _ in range(len(data["chrValues"][chrom]["segments"]))]
        else:
            data["chrValues"][chrom][property] = [0 for _ in range(len(data["chrValues"][chrom]["segments"]))]

        with open(ensemblFn, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in spamreader:
                header = row
                break

            if "Chromosome/scaffold name" in header:
                chromosomeInfoIndex = header.index("Chromosome/scaffold name")
            else:
                chromosomeInfoIndex = header.index("Chromosome/scaffold Name")
            
            aInfoIndex = -1 #index for value with start segment locus
            for i, name in enumerate(header):
                if "tart (bp)" in name:
                    aInfoIndex=i
                    break
            bInfoIndex = -1 #index for value with end segment locus
            for i, name in enumerate(header):
                if "nd (bp)" in name:
                    bInfoIndex=i
                    break
            
            featureNameIndex = header.index(featureName)
            #Found indeces of interesting columbs in the header

            if chromosomeInfoIndex==-1 or aInfoIndex==-1 or bInfoIndex==-1 or featureNameIndex==-1 :
                raise Exception("field not found in header of ensembl file")
            
            #start processing file row by row
            for row in spamreader:
                ch=row[chromosomeInfoIndex]
                ch="chr"+ch
                if ch!=chrom:
                    continue

                A=row[aInfoIndex]
                B=row[bInfoIndex]
                if boolable:
                    #Feature had a name.. then it gets an int
                    feature=row[featureNameIndex]
                    if feature not in data["featureToInt"]:
                        data["featureToInt"][feature] = len(data["featureList"])
                        data["featureList"].append(feature)
                    feature = data["featureToInt"][feature]
                    possibleValues.add(feature)

                            
                A,B = int(A),int(B)
                if B<A: A,B=B,A
                
                overlappingSegments = tree[A:B]
                
                if len(overlappingSegments)>0:
                    iii=9
                    overlappingSegments=list(overlappingSegments)
                    for ovseg in overlappingSegments:
                        lenFeature = B-A
                        lenSegment = ovseg[1]-ovseg[0]
                        overlap = getOverlap(A,B,ovseg[0],ovseg[1])
                        if overlap<0.2*(min(lenFeature,lenSegment)):
                            continue #This overlap is too small
                        if boolable:
                            if feature not in data["chrValues"][chrom][property][ovseg[2]]:
                                data["chrValues"][chrom][property][ovseg[2]].append(feature)
                        else:
                            data["chrValues"][chrom][property][ovseg[2]]+=1
                        
        return data
    def processAllEnsemblFiles(data, ensFns, edir=f"{dataDir}/ensembl/"):
        #Data is a json object
        #ensFns is calculated before, list of file names in edir
        data["featureToInt"] = {}
        data["featureList"] = []
        chrom = data["chrNames"][0]
        
        print(chrom)
        tree = IntervalTree()
        for i, segment in enumerate(data["chrValues"][chrom]["segments"]):
            if len(segment)>0: tree[segment[0]:segment[1]]=i
        
        for ensFileName in ensFns:
            print(ensFileName)
            fullFn = edir+ensFileName
            property = ensToName[ensFileName]
            featureName = ensFnToProperty[ensFileName]
            possibleValues = set()
            boolable = ensFileName in boolableFNs
            processFeatureFile( ensemblFn=fullFn,
                                data=data,
                                tree=tree,
                                property=property,
                                featureName=featureName,
                                possibleValues=possibleValues, 
                                chrom=chrom, boolable=boolable
                            )
            iii=9
        return data

    return processAllEnsemblFiles(data, ensFns, edir=edir) 

#---------------------------------------------------------------------
#Utility function that takes a file, a function that returns [ch,A,B,property] 
#Also takes a list of segemnts.
#Returns another list with corresponding porperties for each segment.
#Afterwards, user may do whatever they want. Either save as a fina list or call this for each tissue type and aggregate the results
#Another param - countabl - bool. If true, just count the rows
#And the chromosome
def getPropertiesForListOfSegments(fileName, segmentList, F, ch, countable=False):
    if fileName=="": 
        #print ("File name empty")
        return [[] for _ in range(len(segmentList))]
    tree = IntervalTree()
    for i, segment in enumerate(segmentList):
        if len(segment)>0: tree[segment[0]:segment[1]]=i
    
    def getOverlap(min1, max1, min2, max2):
        #given 2 segments, return their overlap
        return max(0, min(max1, max2) - max(min1, min2))
    
    ##read the encode datafile
    #Its filename is in ffnn

    
    df = pd.read_csv(fileName, compression='gzip', header=None, sep='\t', dtype=str)
    data_as_list = df.values.tolist()

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
            iii=9
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
    #print(tissue, " searching file name in ", tissueToIDMappingFn)
    df = pd.read_csv(dir+tissueToIDMappingFn, usecols=["Code", "File_ID"], sep='\s*,\s*', engine='python')
    # print(df.iloc[3].Code)
    # print(df.head())
    iii=0
    dataRow = df[df['Code'] == tissue]
    if dataRow.empty:
        print(f"No data file found with {tissue} tissue")
        return []
    else:
        ffnn = dataRow["File_ID"].values[0]
    #ffnn="ENCFF733BFV" #################### BECAUSE MAPPING FILES ARE WRONG!
    # print(ffnn, type(ffnn)) #If empty, then return "", processed later, all props will be empty for the tissue
    if type(ffnn)==float: return ""
    
    
    ffnn = dir+ffnn+".bed.gz"
    return ffnn

def encodeGetter(row):
    #col_names = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "featureType", "classification"]
    return ([row[0], int(row[1]), int(row[2]), row[9].split(',')]) #TODO SERIOUSLY?

#Adds information from 
def addEncode(data, tfn=f"{encodeTissueMappingFn}", dir=f"{dataDir}/ENCODE/"):
    ch=data["chrNames"][0]
    data["chrValues"][ch]["encodeStates"]=[dict() for _ in range(len(data["chrValues"][ch]["segments"]))]
    if "tissueSpecificProperties" not in data: data["tissueSpecificProperties"]=[]
    if "encodeStates" not in data: data["tissueSpecificProperties"].append("encodeStates")
    possibleValues = set()
    U=data
    for tis in list(U['tissueData']["tissueNames"].keys()):
        if tis=="AO":
            iii=9
        bit = U['tissueData']['tissueBits'][tis]
        links = U["chrValues"][ch]["links"]
        segments = U["chrValues"][ch]["segments"]
        tisLinks = [link for link in links if (link[2]&bit)==bit]
        tisSegments = [link[0] for link in tisLinks]+[link[1] for link in tisLinks]
        tisSegments = set(tisSegments)
        segs = [segment if i in tisSegments else [] for i,segment in enumerate(segments)]

        RR = getPropertiesForListOfSegments(getFileName(tis,tfn,dir), segs, F=encodeGetter, ch=ch, countable=False)
        for i,segmentSetOfProps in enumerate(RR):
            if len(segmentSetOfProps)==0: continue
            data["chrValues"][ch]["encodeStates"][i][tis]=list(segmentSetOfProps)
            for el in segmentSetOfProps: possibleValues.add(el)
    data["possibleEncodes"] = list(possibleValues)
    return data

#-----------------------------------------------
#addGeneGTExData, taken and adapted from addGeneDistance.ipynb
def addGeneGTExData(U):
    import csv, json
    GTExfn = f"{dataDir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
    GRCh38fn = f"{dataDir}/ensembl_genes_GRCh38.txt"
    GTExTissueTranslationFn = f"{dataDir}/{GTExTissueMappingFn}"
    #read GTEx file to G
    with open(GTExfn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        G = list(spamreader)
        G = G[2:] #remove first 2 lines
        headerG = G.pop(0) #remove the header

    def parse_gene_info(gene_info_str):
        # Split the string by commas
        parts = gene_info_str.split(',')

        # Assign each part to a variable
        gene_id = parts[0]
        transcript_id = parts[1]
        chromosome = parts[2]
        [start_position, end_position] = sorted([int(parts[3]), int(parts[4])])
        strand = int(parts[5])
        gene_name = parts[6]
        gene_type = parts[7]

        # Return a dictionary with the parsed data
        return {
            'gene_id': gene_id,
            'transcript_id': transcript_id,
            'chromosome': "chr"+chromosome,
            'start': start_position,
            'end': end_position,
            'strand': strand,
            'gene_name': gene_name,
            'gene_type': gene_type
        }


    #read the ensembl_genes_GRCh38.txt
    with open(GRCh38fn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        T = list(spamreader)
        headerT = T.pop(0) #remove the header
        T = [parse_gene_info(el[0]) for el in T]

    #make translation easier. Structure TT[ch][geneName] gets the [start, end]
    TT = {}
    for dic in T:
        ch = dic["chromosome"]
        #if len(ch)>2: continue #get rid of weird chromosome names
        ch="idk"
        if ch not in TT: TT[ch]={}
        if dic["gene_id"] in TT[ch]:
            raise Exception("duplicate gene_id found")
        TT[ch][dic["gene_id"]] = sorted([dic["start"], dic["end"]])+[dic["chromosome"]]

    #read file with tissue to column name
    with open(GTExTissueTranslationFn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        sprl = list(spamreader)
        tisToTitleRows = {row[0].strip():row[1] for row in sprl}
        titleRowsToTis = {row[1]:row[0].strip() for row in sprl}

    #create a list with objects about each row of G
    #{rowInd:0, ensg: "", start:1000, end:1100, chr:"chr7", tissues:set(AD, AO)}
    Ginfo = []
    for i in range(len(G)):
        ob = {}
        ob["ensg"] = G[i][0].split('.', 1)[0] #remove .5 or .6 or... after id
        ob["ensgName"] = G[i][1]
        ensg = ob["ensg"]
        if ensg not in TT["idk"]: continue #this row will be ignored
        [ob["start"],ob["end"],ob["chr"]] = TT["idk"][ob["ensg"]]
        ob["tissues"] = set()
        ob["tissueExpressions"] = {}
        for j in range(len(headerG)):
            if j<2: continue
            #if float(G[i][j])<=0: continue #change - dont ignore zero expressions
            #only non zero expressions used to remain, not they remain too
            curTissue = headerG[j] #e.g. Adipose - Subcutaneous
            if curTissue not in titleRowsToTis:
                continue
            curTis = titleRowsToTis[curTissue] #e.g. AD2
            ob["tissues"].add(curTis)
            ob["tissueExpressions"][curTis] = float(G[i][j])
        if len(ob["tissues"])==0: continue #not a single of our tissues has >0 expression


        Ginfo.append(ob)
    
    def getMinDistancesForAllSegments(segments, fantoms):
        i=0
        j=0 # currentHypothesisForClosestFantom
        def skelas(a,b):
            return a[0]<=b[1] and b[0]<=a[1]  #!!!
        def isOnTheLeft(star, makonis): #star is a fantom segment, makonis is a segment. 
            return star[1]<makonis[0]
        fantoms.append([10e14, 10e15+10,-1]) #append a dummy fantom that is definetly not the closest one to any segment
        res = [10e10 for _ in range(len(segments))] #for each segment and for each tissue
        fantomIndeces = [[] for _ in range(len(segments))] #for each segment, find the name of the fantom, use 2nd column of the data file
        while (i<len(segments)) :
            makonis = segments[i]
            star = fantoms[j]

            if skelas(makonis, star):
                
                if res[i]>0: #ja tur jau bija nenulle
                    fantomIndeces[i]=[] #ja pirms tam bija kaut kas aizpildīts un tgd ieraugām kaut ko, kas šķeļas, to apnuļļo
                res[i]=0
                fantomIndeces[i].append(star[2])

                J=j
                STAR=star
                while(J<len(fantoms)-1):
                    J+=1
                    STAR=fantoms[J]
                    if skelas(makonis, STAR):
                        fantomIndeces[i].append(STAR[2])
                    else: 
                        break
                i+=1
                continue
            
            if (isOnTheLeft(star, makonis)):
                res[i]=min(res[i],makonis[0]-star[1])
                fantomIndeces[i]=[star[2]]
                j+=1
                continue
            dist = star[0]-makonis[1]
            
            if j>0 and makonis[0]-fantoms[j-1][1] < dist :
                #previous fantom was closer
                dist=makonis[0]-fantoms[j-1][1]
                fantomIndeces[i]=[fantoms[j-1][2]]
            else:
                fantomIndeces[i]=[fantoms[j][2]]
            res[i] = dist

            i+=1
        fantoms.pop(-1) #remove the dummy element from the back. ~10:10
        return [res,fantomIndeces]



    def getGeneExpressionValue(rowId, tissue, Src):
        if rowId<0:
            raise Exception(f"rowId can not be negative {rowId}")
        if tissue in Src[rowId]["tissues"]:
            return Src[rowId]["tissueExpressions"][tissue]
        else: return 0

    def getGeneName(rowId, Src):
        if rowId<0:
            raise Exception(f"rowId can not be negative {rowId}")
        return Src[rowId]["ensgName"] 

    if True:
        iii=0
        #Finding bitmaps for segments
        ch = U["chrNames"][0]
        # allTriangles = U["triangles"] 
        allSegments = U["chrValues"][ch]["segments"]
        allLinks = U["chrValues"][ch]["links"]
        
        segmentBitmaps = [0 for _ in range(len(allSegments))]
        for [a,b,bit] in allLinks:
            segmentBitmaps[a]|=bit
            segmentBitmaps[b]|=bit
        
        curChrDataRows=[]           
        curChrDataRows = [el for el in Ginfo if el["chr"]==ch]
        curChrDataRows = sorted(curChrDataRows, key=lambda x: (x["start"], x["end"]))

        fantomResults = [{} for _ in range(len(allSegments))]
        fantomValues = [{} for _ in range(len(allSegments))]
        fantomNames = [{} for _ in range(len(allSegments))]

        for tissue in U["tissueData"]["tissueBits"].keys():
            tisBit = U["tissueData"]["tissueBits"][tissue]
            tissueSegments = [allSegments[i] for i in range(len(allSegments)) if ((segmentBitmaps[i]&tisBit)==tisBit)]
            tissueSegmentInds = [i for i in range(len(allSegments)) if ((segmentBitmaps[i]&tisBit)==tisBit)]
        
            fantomSegments = \
                [[curChrDataRows[i]["start"],curChrDataRows[i]["end"],i] for i in range(len(curChrDataRows)) if tissue in curChrDataRows[i]["tissues"]]#

            [distToF, bestGeneIndsInCurChrDataRows] = getMinDistancesForAllSegments(tissueSegments, fantomSegments)
            iii=0

            for k,dist in enumerate(distToF):
                #k is index, dist is the distance from segment to the closest fantom
                #store results to fantomResults list of dicts
                segIndInAllSegments = tissueSegmentInds[k] #distToF only has data on segments that are in the current tissue. Tranlating the index to index in allSegments 
                fantomResults[segIndInAllSegments][tissue] = distToF[k]

                fantomValues[segIndInAllSegments][tissue] = max([getGeneExpressionValue(ind, tissue, curChrDataRows) for ind in bestGeneIndsInCurChrDataRows[k]])  #bestExprs[k] 
                fantomNames[segIndInAllSegments][tissue] = [getGeneName(ind,curChrDataRows) for ind in bestGeneIndsInCurChrDataRows[k]]#
                
            iii=0 #end processing one tissue

        iii=0 #end processing all tissues
        U["chrValues"][ch]["closestGTExGenes"] = fantomResults
        U["chrValues"][ch]["segmentGTExGenes"] = fantomValues
        U["chrValues"][ch]["GTExGeneNames"] = fantomNames
        print("addGeneGTExData done")
        return U
        #######################################################################################
        # ffnn = udir+f"wGenes15-{ch}.json"
        # with open(ffnn, 'w', encoding ='utf8') as json_file: 
        #     json.dump(U, json_file, ensure_ascii = False) 
        # print(ffnn)
        
def encodeRNAGetter(row):
    if (float(row[8])<1): return []
    return [row[0],int(row[1]), int(row[2]), 1]

def addEncodeRNA(data, tfn=ENCODE_RNA_tissue_mappingFn, dir=f"{dataDir}/ENCODE/"):
    ch=data["chrNames"][0]
    data["chrValues"][ch]["encodeRNAStates"]=[dict() for _ in range(len(data["chrValues"][ch]["segments"]))]
    if "tissueSpecificProperties" not in data: data["tissueSpecificProperties"]=[]
    if "encodeRNAStates" not in data: data["tissueSpecificProperties"].append("encodeRNAStates")
    
    U=data
    for tis in list(U['tissueData']["tissueNames"].keys()):
        if tis=="AO":
            iii=9
        bit = U['tissueData']['tissueBits'][tis]
        links = U["chrValues"][ch]["links"]
        segments = U["chrValues"][ch]["segments"]
        tisLinks = [link for link in links if (link[2]&bit)==bit]
        tisSegments = [link[0] for link in tisLinks]+[link[1] for link in tisLinks]
        tisSegments = set(tisSegments)
        segs = [segment if i in tisSegments else [] for i,segment in enumerate(segments)]

        RR = getPropertiesForListOfSegments(getFileName(tis,tfn,dir), segs, F=encodeRNAGetter, ch=ch, countable=True)
        for i,propCount in enumerate(RR):
            if propCount==0 or type(propCount)==type([]): continue
            data["chrValues"][ch]["encodeRNAStates"][i][tis]=propCount
    return data
    

import networkx as nx
from networkx.exception import PowerIterationFailedConvergence


def addCentralities(U):
    #Finding bitmaps for segments
    ch = U["chrNames"][0]
    allSegments = U["chrValues"][ch]["segments"]
    U["chrValues"][ch]["ECentrality"]=[{} for _ in range(len(allSegments))]
    U["chrValues"][ch]["CCentrality"]=[{} for _ in range(len(allSegments))]
    
    for tissue in U["tissueData"]["tissueBits"].keys():
        tisBit = U["tissueData"]["tissueBits"][tissue]
        tissueLinks = [(a,b) for [a,b,bit] in U["chrValues"][ch]["links"] if ((bit&tisBit)==tisBit)]
        tissueLinks = sorted(list(set(tissueLinks)))
        iii=0

        if CALCULATE_CENTRALITIES:
            G = nx.Graph()
            G.add_edges_from(tissueLinks)
            # Calculate eigenvector centrality
            try:
                centrality = nx.eigenvector_centrality(G, max_iter=2000)
            except PowerIterationFailedConvergence as e:
                centrality = {node: 1 for node in G.nodes()}
            # Save the centrality of each node
            for node, centrality_value in centrality.items():
                U["chrValues"][ch]["ECentrality"][node][tissue] = centrality_value
                #print(f"Node {node}: {centrality_value}")
            centrality = nx.closeness_centrality(G)
            for node, centrality_value in centrality.items():
                U["chrValues"][ch]["CCentrality"][node][tissue] = centrality_value
        else:
            for node in range(len(U["chrValues"][ch]["ECentrality"])):
                U["chrValues"][ch]["ECentrality"][node][tissue] = -1
                U["chrValues"][ch]["CCentrality"][node][tissue] = -1

def makeAdjAndBitsOfLink(links):
    # creates an adj structure for all links and creates a dictionary bitsOfLink
    # bitsOfLink[(A,B)] == the bitmap of the link
    adj = dict()
    bitsOfLink = dict()
    for link in links:
        [A, B, bitmap, *pv] = link
        if (B<=A): 
            #print ("Warning: got an unsorted Link in makeAdj method!")
            #print(link)
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
    # uses links == self.links
    # links argumentthat is used here has list of links in form [A, B, bitmap] for one chromosome
    
    [adj, bitsOfLinks] = makeAdjAndBitsOfLink(links) 
    #adj[A] == set of all other vertices that are adjacent to A
    #bitsOfLink[(A,B)] == the bitmap of the link (A,B)
    cliques = set() #{(A,B,C), (), ()}
    
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
def mkVizFile(data):
    ch = data["chrNames"][0]
    cc = getC3(data["chrValues"][ch]["links"])
    #Cik ir virsotnes un cik ir linki in cc
    links = set()
    segs = set()
    for (A,B,C,bits) in cc:
        for tup in [(A,B), (B,C), (A,C)]:
            links.add(tup)
        for seg in [A,B,C]:
            segs.add(seg)
    print("C3 calculated")
    newSegs = sorted(list(segs)) #Saraksts no veciem indeksiem, kas paliek, jo ir trissturos
    newSegIndex = {newSegs[i]: i for i in range(len(newSegs))} #Oriģinālajam indeksam piekārto jauno indeksu 
    oldSegIndex = {i: newSegs[i] for i in range(len(newSegs))} #Jaunajam indeksam piekārto oriģinālo
    newData = deepcopy(data)
    newData["segments"] = [data["chrValues"][ch]["segments"][oldSegIndex[i]] for i in range(len(newSegs))]
    newTriangles = [(newSegIndex[A], newSegIndex[B], newSegIndex[C], bits) for (A,B,C,bits) in cc]
    newData["triangles"] = newTriangles


    L = []

    possibleValsFor={}
    for tag in data["boolable properties"]:
        possibleValsFor[tag]=set()
        for lst in data["chrValues"][ch][tag]:
            for el in lst:
                possibleValsFor[tag].add(el)

    iii=9
    print(possibleValsFor.keys())
    nickname = {
        "miRNA_Target_Regions": "miRNA targ",
        "Other_Regulatory_Regions": "Regul. reg",
        "Regulatory_Features": "Reg. Feat",
        'Somatic_Short_Variants': 'Som. sh. v',
        'Somatic_Structural_Variants': 'Som. str. v',
        'Structural_Variants': 'Str. v',
    }
    iii=9

    for i in range(len(newData["segments"])):
        header = []
        row = [i, newData["segments"][i][0], newData["segments"][i][1]]
        header=["ID", "A", "B"]
        for tag in data["boolable properties"]:
            for propID in possibleValsFor[tag]:
                prop = nickname[tag]+":"+data["featureList"][propID]
                header.append(prop)
                row.append(1 if propID in data["chrValues"][ch][tag][oldSegIndex[i]] else 0)
        for tag in data["countable properties"]:
            prop=nickname[tag]+":"+"count"
            header.append(prop)
            row.append(data["chrValues"][ch][tag][oldSegIndex[i]])

        L.append(row)
    ii=9



    newData["segmentEnsemblHeader"] = header
    newData["segmentEnsembl"] = L

    iii=9
    newData["chrValues"][ch]["segmentStates"] = [data["chrValues"][ch]["segmentStates"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["encodeStates"] = [data["chrValues"][ch]["encodeStates"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["segmentGenes"] = [data["chrValues"][ch]["segmentGenes"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["encodeRNAStates"] = [data["chrValues"][ch]["encodeRNAStates"][oldSegIndex[i]] for i in range(len(newSegs))]
    #HERE TODO?
    #Updating closestFantoms
    newData["chrValues"][ch]["closestFantoms"] = [data["chrValues"][ch]["closestFantoms"][oldSegIndex[i]] for i in range(len(newSegs))]
    #Updating segmentFantoms
    newData["chrValues"][ch]["segmentFantoms"] = [data["chrValues"][ch]["segmentFantoms"][oldSegIndex[i]] for i in range(len(newSegs))]
    #Updating fantomNames
    newData["chrValues"][ch]["fantomNames"] = [data["chrValues"][ch]["fantomNames"][oldSegIndex[i]] for i in range(len(newSegs))]
    
    newData["chrValues"][ch]["ECentrality"] = [data["chrValues"][ch]["ECentrality"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["CCentrality"] = [data["chrValues"][ch]["CCentrality"][oldSegIndex[i]] for i in range(len(newSegs))]
    
    newData["chrValues"][ch]["closestGTExGenes"] = [data["chrValues"][ch]["closestGTExGenes"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["segmentGTExGenes"] = [data["chrValues"][ch]["segmentGTExGenes"][oldSegIndex[i]] for i in range(len(newSegs))]
    newData["chrValues"][ch]["GTExGeneNames"] = [data["chrValues"][ch]["GTExGeneNames"][oldSegIndex[i]] for i in range(len(newSegs))]




    for tag in data["boolable properties"]+data["countable properties"]+["segments"]:
        if tag in newData["chrValues"][ch]:
            del newData["chrValues"][ch][tag]
    

    
    newData["VersionInfo"] = {"Generated using": "processEncodeEnsembl.py",
                              "inputFile": f"{dataDir}/wfantoms27-{ch}.json",
                              "inputFileMadeBy": "processEncodeEnsembl.py",
                              "purpose": "File ready for visualization",
                              "version": "05-02-24-1",
                              }
    print("Start dumping file")
    fffnn=f"{udir}/Visualizable-{DS}-{ch}-cliques-genes-encode-GTEx-Fantoms-centr-viz0205.json"
    with open(fffnn, 'w') as file:
        json.dump(newData, file) 
    print("Done saving file", fffnn)


def processUniversal(ufn):
    with open(ufn, 'r') as file:
        dataU = json.load(file)
    dataU["tissueSpecificProperties"] = ["segmentStates"]
    data = translateGeneFile(dataU) #It is now translated to hg38
    
    #get rid of invalid segments and update link indeces
    #Also update segmentGenes and segmentStates lists
    data=validifyUniversal(data)
    print("validifying done")

    data=addEncodeRNA(data)
    print("EncodeRNA added")

    data = addEncode(data)
    print("Encode added")
    

    data = processAllEnsembl(data)
    print("Ensembl processed")

    data = addGeneGTExData(data)

    #calc centrality
    addCentralities(data)
    
    mkVizFile(data)
    
    return data

def prepareSummary(S, fn):
    ch = S["chrNames"][0]
    chData = S["chrValues"][ch]
    C3 = getC3(chData["links"])
    bigRez = [] #for all tissues

    for tis in S["tissueData"]["tissueIDs"]:
        #tis is e.g. 'PA'
        tisBit = S["tissueData"]["tissueBits"][tis]
        tisLinks = [[A,B] for [A,B,bit] in chData["links"] if (bit&tisBit)==tisBit]
        tisSegments = set([A for [A,B] in tisLinks]+[B for [A,B] in tisLinks])
        
        tisCliques = [(A,B,C,bit) for (A,B,C,bit) in C3 if (bit&tisBit)==tisBit]
        tisCliqueSegments = [A for (A,B,C,bit) in tisCliques]+[B for (A,B,C,bit) in tisCliques]+[C for (A,B,C,bit) in tisCliques]
        tisCliqueSegments = set(tisCliqueSegments)

        if len(tisCliqueSegments.difference(tisSegments))>0:
            raise Exception("clique segment set is not a subset of all tissue segment set")

        REZ = []
        MAXALLOWEDDISTANCE=10e15
        for tisSegment in tisSegments:
            record = {
                        "chr": ch,
                        "tissue": tis,
                        "tissueName": S["tissueData"]["tissueNames"][tis],
                        "segmentID": tisSegment,
                        "start": chData["segments"][tisSegment][0],
                        "end": chData["segments"][tisSegment][1],
                    }
            iii=0
            try:
                closestFantomInfo = {
                    "names": chData["fantomNames"][tisSegment][tis],
                    "distance": chData["closestFantoms"][tisSegment][tis],
                    "maxExpression":chData["segmentFantoms"][tisSegment][tis], 
                }
                closestGTExInfo = {
                    "names": chData["GTExGeneNames"][tisSegment][tis],
                    "distance": chData["closestGTExGenes"][tisSegment][tis],
                    "maxExpression":chData["segmentGTExGenes"][tisSegment][tis], 
                }
            except KeyError:
                raise Exception("key not found in dict, should have been found")
            
            record["isInClique"] = (tisSegment in tisCliqueSegments)

            record["fantomExpression"]="NA" if (closestFantomInfo["distance"]>MAXALLOWEDDISTANCE) else closestFantomInfo["maxExpression"]
            record["fantomDistance"] = closestFantomInfo["distance"]
            record["fantomCount"] = len(closestFantomInfo["names"])

            record["GTExExpression"]="NA" if (closestGTExInfo["distance"]>MAXALLOWEDDISTANCE) else closestGTExInfo["maxExpression"]
            record["GTExDistance"] = closestGTExInfo["distance"]
            record["GTExCount"] = len(closestGTExInfo["names"])

            record["RNAPolCount"] = 0 if tis not in chData["encodeRNAStates"][tisSegment] else chData["encodeRNAStates"][tisSegment][tis]
            iii=0

            REZ.append(record)
            iii=0
        iii=0
        bigRez.extend(REZ)

    iii=9

    cols = list(REZ[0].keys())

    table = []
    table.append(cols) #add header
    for record in bigRez:
        vals = [record[colName] for colName in cols]
        table.append(vals)

    with open(fn, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerows(table)


#Call
for ch in CHRS:
    print("<<<<<<<<<<<<<<<<<<<<<<<", ch)
    data1 = processUniversal(f"{udir}/rez/wfantoms05-{ch}.json")
    with open(f"{udir}/rez/universal-{DS}-with-extras-{ch}-hg39-0205.json", 'w') as file:
        json.dump(data1, file) 
    ffknn = f"{udir}/rez/overview-Fantoms-GTEx-Segments-{DS}-{ch}-v0205.csv"
    prepareSummary(data1, ffknn)

