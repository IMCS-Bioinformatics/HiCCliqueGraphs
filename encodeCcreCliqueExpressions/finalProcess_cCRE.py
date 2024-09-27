import json, csv
import itertools
import gc
import os
import csv

#######################################################################################
## This module creates a .csv file with information about expressions of cCRE cliques
## Input: JSON file created by prepareForcCRE.py
## Output: csv file
#######################################################################################

#Change below to choose dataset !!!!!!!!!!!!!!!!!!!!!!!!!!!
DS = "tisPCHiC-mini" #Tissue pcHi-C dataset with chr5 and chr10 only, included in repository
# DS = "tisPCHiC" #Tissue pcHi-C
# DS = "3DIV" #Tissue Hi-C

# Directory with the current .py module
currentDir = os.path.dirname(os.path.abspath(__file__))


if DS=="tisPCHiC":
    Ufn = "encodeCcreCliqueExpressions/results/data-pvalue-0.7-fin-digraph-GRCh38-wEncodeBits-wExpressions-wC3.json"
    rezFn = "encodeCcreCliqueExpressions/results/stats-encodeCCRE-tisPCHiC.csv"
elif DS=="tisPCHiC-mini":
    Ufn = "encodeCcreCliqueExpressions/results/data-pvalue-0.7-fin-digraph-GRCh38-mini-wEncodeBits-wExpressions-wC3.json"
    rezFn = "encodeCcreCliqueExpressions/results/stats-encodeCCRE-tisPCHiC-mini.csv"
elif DS=="3DIV":
    Ufn = "encodeCcreCliqueExpressions/results/data-pvalue-10-fin-digraph-GRCh38-wEncodeBits-wExpressions-wC3.json"
    rezFn = "encodeCcreCliqueExpressions/results/stats-encodeCCRE-3DIV.csv"

headerMappingFN = f'{currentDir}/header_mapping.json'

#Read data file
with open(Ufn, 'r') as f:
    U = json.load(f)

# Utility function to check if clique atisfies a condition
def cliqueHasOneSameProperty(clique, encodeProperty, ch, tis):
    propertyBit = U["encodeBits"][encodeProperty]
    encodesForSegments = U["chrValues"][ch]["encodeStates"]
    if len(   [v for v in clique if (tis in encodesForSegments[v]) and ((encodesForSegments[v][tis]&propertyBit)==propertyBit)  ])==3:
        return True
    return False

# Generate all possible combinations of length 3 with repetition allowed
letters = ['P', 'E', 'C', 'D']
combinationsPECD = list(itertools.combinations_with_replacement(letters, 3))

getFullName = {'P':"promoter",
               'E': "distal enhancer",
               'C': "CTCF",
               'D':"DNase-only"}

#Utility function to make a str from a list
def getTripletStr(triplet):
    s= ""
    for el in triplet:
        s+=str(el)
    return s

#Utility function to check if clique has these exact 3 cCRE properties in the nodes
def cliqueSatisfiesEncodeTriplet(clique, ch, tis, triplet):
    #Triplet - tuple of three elements, each is a letter
    triplets = set(list(itertools.permutations(triplet)))
    for tr in triplets:
        tmp = [w for i,w in enumerate(tr) if (U["chrValues"][ch]["encodeStates"][clique[i]].get(tis, 0) & U["encodeBits"][getFullName[w]])>0]
        if len(tmp)==3:
            return True
    return False

import statistics, numpy as np
statisticsProps = ["mean","25%", "50%", "75%", "cardinality"]
#Utility function to calculate expression statistics for a set of expression values
def getBasicStatistics(data):
    if len(data)==0:
        return {k:-1 for k in statisticsProps}
    mean = statistics.mean(data)
    # Calculate quartiles
    quartiles = np.percentile(data, [25, 50, 75])
    return {
        "mean": mean,
        "25%": quartiles[0],
        "50%": quartiles[1],  # This is the median
        "75%": quartiles[2],
        "cardinality": len(data) #number of genes or FANTOM5 profiles, i.e. how many expression values were available
    }

#Calculates expression statistics (both genes and fantom5) for a set of segments
def getExpressionsForSegments(setOfSegments, ch, tis):
    #returns statistics about expressions, as a dict
    props = statisticsProps
    res = {
        "GTEx":{k:0 for k in props},
        "FANTOM5":{k:0 for k in props},
    }
    #Collect all gene IDs

    ###Process GTEx
    if tis not in U["dataAviabilityForTissues"]["GTEx"]:
        res["GTEx"] = {k:-1 for k in props}
    else:
        setOfActiveGenes = set()
        for segID in setOfSegments:
            setOfActiveGenes.update(U["chrValues"][ch]["segmentGenes"][segID])
        geneExpressions = [ U["chrValues"][ch]["chromosomeGenes"][geneID]["GTEx"].get(tis, -1) for geneID in setOfActiveGenes]
        geneExpressions = [ge for ge in geneExpressions if ge!=-1]
        stats = getBasicStatistics(geneExpressions)
        for k in props:
            res["GTEx"][k] = stats[k]
    
    ###Process FANTOM5
    if tis not in U["dataAviabilityForTissues"]["fantoms"]:
        res["FANTOM5"] = {k:-1 for k in props}
    else:
        setOfActiveFantoms = set()
        for segID in setOfSegments:
            setOfActiveFantoms.update(U["chrValues"][ch]["segmentFantoms"][segID])
        fantomExpressions = [ U["chrValues"][ch]["chromosomeFantoms"][fantomID]["expressionInTissues"].get(tis, -1) for fantomID in setOfActiveFantoms]
        fantomExpressions = [ge for ge in fantomExpressions if ge!=-1]
        stats = getBasicStatistics(fantomExpressions)
        for k in props:
            res["FANTOM5"][k] = stats[k]
    return res

# Modifies U object, adding data for cliques of size up to maxR, for one tissue
# Data is saved in  U["chrValues"][ch][f"cliques{RR}"] and U["chrValues"][ch][f"cliques{RR}_minus1"]
# Calculated data will be replaced every time this method is called, to save memory 
# Almost cliques are defined as cliques without one missing edge
def addCliquesAndAlmostCliques(U,ch, tis, maxR=4,):
    cliques = U["chrValues"][ch]["cliques"] #Cliques of size R are calculated, using cliques of size R-1
    for R in range(4,maxR+1):
        cliquesX = addC4ToU2(U,ch, cliques, R, tis) #returns a list of full cliques of size R
        cliques=cliquesX

#Utility function that checks the results and raises exception if any clique that was found is not a valid clique
def checkCliques(U, ch, cliques, checkFullClique=True):
    #clique - list or tuple where last element is the bitmap of tissues
    edgeNodeInds = list(itertools.combinations(list(range(0,len(cliques[0])-1)), 2))
    c=0

    for tis, tisBit in U["tissueBits"].items():
        #Checks cliques of each tissue type seperately
        tissueCliques=  [clique for clique in cliques if ((clique[-1]&tisBit)==tisBit)]
        tissueLinks=  [link for link in U["chrValues"][ch]["links"] if ((link[2]&tisBit)==tisBit)]
        adjTis = {}
        for link in tissueLinks:
            if link[0] not in adjTis:
                adjTis[link[0]]=set()
            if link[1] not in adjTis:
                adjTis[link[1]]=set()
            adjTis[link[0]].add(link[1])
            adjTis[link[1]].add(link[0])
        #Check all cliques

        for clique in tissueCliques:
            c = clique[:-1] #remove bit
            existingLinkCountOfClique = 0
            for edge in edgeNodeInds: #(0,1), (0,2), ..., (2,3) is edge
                if (clique[edge[1]] in adjTis[clique[edge[0]]]) and (clique[edge[0]] in adjTis[clique[edge[1]]]):
                    existingLinkCountOfClique+=1
            if checkFullClique:
                if existingLinkCountOfClique==(len(c)*(len(c)-1)//2):
                    pass
                else:
                    raise Exception("Clique identified as full clique does not have all links")
            else:
                #Checking almost clique
                if existingLinkCountOfClique>=(len(c)*(len(c)-1)//2)-1:
                    pass
                else:
                    raise Exception("Clique identified as almost clique does not have all links")
        if c==0: c=[]
        print(f"\t\t{tis} OK, had {len(tissueCliques)} cliques{'_minus1' if (not checkFullClique) else ''} of size {len(c)}")
    print(f"{ch} OK")

#caliucultes cliques of size RR, for chromosome ch and tissue tis. Uses cliques of RR-1 in parameter cliques
def addC4ToU2(U, ch, cliques, RR, tis):
    fullCliquesChr = dict()
    almostCliquesChr = dict()
    tisBit = U["tissueBits"][tis]
    print(f"start calculating cliques for {tis}")
    tissueCliques =[clique[:-1] for clique in cliques if ((clique[-1]&tisBit)==tisBit)]
    tissueLinks = [link for link in U["chrValues"][ch]["links"] if ((link[2]&tisBit)==tisBit)]
    #make adj
    adjTis = dict()
    for link in tissueLinks:
        if link[0] not in adjTis:
            adjTis[link[0]]=set()
        if link[1] not in adjTis:
            adjTis[link[1]]=set()
        adjTis[link[0]].add(link[1])
        adjTis[link[1]].add(link[0])
    
    #Find full cliques of the next size
    #RR will be 4 when I am looking for cliques of size 4
    tissueFullCliques = set()
    for clique in tissueCliques:
        neighborIntersection = adjTis[clique[0]]
        for node in clique[1:]:
            neighborIntersection=neighborIntersection.intersection(adjTis[node])
        #each node in neighborIntersection forms a clique of next size
        for d in neighborIntersection:
            tissueFullCliques.add(tuple(sorted(clique+[d])))
    
    for clique in tissueFullCliques:
        fullCliquesChr[clique]=(fullCliquesChr.get(clique,0)|tisBit)

    #Find all almost cliques
    #The invariant - each almost-clique has one full clique of size RR-1
    #I loop over cliques of size RR-1 again and find nodes that are neighbors to exactly all but one node
    tissueAlmostCliques = set()
    for clique in tissueCliques:
        neighborCounts = {}
        for node in clique:
            for neig in adjTis[node]:
                neighborCounts[neig]=neighborCounts.get(neig,0)+1
        formsAlmostClique = [node for node, count in neighborCounts.items() if (count==RR-2) and node not in clique] #if count==RR-2 then node is neighbor to all but one node => forms an almost clique
        for d in formsAlmostClique:
            tissueAlmostCliques.add(tuple(sorted(clique+[d])))
    print(f"Calculated: {ch}\t{tis}\tC{RR-1}:{len(tissueCliques)}\tC{RR}:{len(tissueFullCliques)}\tC{RR}_minus1:{len(tissueAlmostCliques)}")
    for clique in tissueAlmostCliques:
        almostCliquesChr[clique]=(almostCliquesChr.get(clique,0)|tisBit)


    fullCliquesChr = [list(tup)+[bit] for tup,bit in fullCliquesChr.items()]
    almostCliquesChr = [list(tup)+[bit] for tup,bit in almostCliquesChr.items()]
    U["chrValues"][ch][f"cliques{RR}_minus1"] = almostCliquesChr
    U["chrValues"][ch][f"cliques{RR}"] = fullCliquesChr #Overwrites the results in U

    return fullCliquesChr



#Utility function to get statistics of GTEx and FANTOM5 expressions, accepts list of cliques as input
def getExpressionsForListOfCliques(listOfCliques, ch, tis):
    #calls getExpressionsForListOfSegments
    segs = set([clique[i] for clique in listOfCliques for i in [0,1,2]])
    return getExpressionsForSegments(segs, ch, tis)




##################################################################

#Main function that calculates everything for a given list of chromsomes
#I.e. for each tissue type, for each chromsome calculates expression statistics for different types of nodes and cliques
def calculatecCRE_stats(chrs):
    L=[] #CSV rows will be stored here. Each row for one chr, one tissue type
    R={} #Aggregates info about expressions, could be used to make plots. Currently not used
    for ch in chrs:
        if ch not in U["chrNames"]:
            continue
        print(ch)

        for tis in U["dataAviabilityForTissues"]["encodeStates"]:
            #Only process tissues that have cCRE properties available

            header = ["chr", "tissue"] #Header will be filled iteratively
            row = [ch, tis] #And row will be filled at the same time

            tisBit = U["tissueBits"][tis] #Tissue information is stored as bitmaps

            ##################### "segmnets", "links"
            tissueNodes = [link[0] for link in U["chrValues"][ch]["links"] if ((link[2]&tisBit)==tisBit)]+[link[1] for link in U["chrValues"][ch]["links"] if ((link[2]&tisBit)==tisBit)]
            tissueNodes = set(tissueNodes)
            tissueLinks = [link[:2] for link in U["chrValues"][ch]["links"] if ((link[2]&tisBit)==tisBit)]
            header+=["segments in graph", "links in graph"]
            row+=[len(tissueNodes), len(tissueLinks)]

            curStats = getExpressionsForSegments(tissueNodes, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} {typ}: for all nodes"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("all nodes", []).append(curStats)
            

            ##################### P segments in graph with promoter
            segmentsWithP = [nodeID for nodeID in tissueNodes if (U["chrValues"][ch]["encodeStates"][nodeID].get(tis, 0) & U["encodeBits"]["promoter"])==U["encodeBits"]["promoter"]]
            segmentsWithP=set(segmentsWithP)
            header+=["P: segments in graph with promoter"]
            row+=[len(segmentsWithP)]
            
            curStats = getExpressionsForSegments(segmentsWithP, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} P {typ}: for all segments with promoters"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("promoter nodes", []).append(curStats)

            segmentsWithE = [nodeID for nodeID in tissueNodes if (U["chrValues"][ch]["encodeStates"][nodeID].get(tis, 0) & U["encodeBits"]["distal enhancer"])==U["encodeBits"]["distal enhancer"]]
            segmentsWithE=set(segmentsWithE)
            header+=["E: segments in graph with distal enhancer"]
            row+=[len(segmentsWithE)]
            
            curStats = getExpressionsForSegments(segmentsWithE, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} E {typ}: for all segments with distal enhancers"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("d.enhancer nodes", []).append(curStats)


            segmentsWithC = [nodeID for nodeID in tissueNodes if (U["chrValues"][ch]["encodeStates"][nodeID].get(tis, 0) & U["encodeBits"]["CTCF"])==U["encodeBits"]["CTCF"]]
            segmentsWithC=set(segmentsWithC)
            header+=["C: segments in graph with CTCF"]
            row+=[len(segmentsWithC)]
            
            curStats = getExpressionsForSegments(segmentsWithC, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C {typ}: for all segments with CTCF"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("CTCF nodes", []).append(curStats)

            segmentsWithD = [nodeID for nodeID in tissueNodes if (U["chrValues"][ch]["encodeStates"][nodeID].get(tis, 0) & U["encodeBits"]["DNase-only"])==U["encodeBits"]["DNase-only"]]
            segmentsWithD=set(segmentsWithD)
            header+=["D: segments in graph with DNase-only"]
            row+=[len(segmentsWithD)]
            
            curStats = getExpressionsForSegments(segmentsWithD, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} D {typ}: for all segments with DNase-only"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("DNase-only nodes", []).append(curStats)
            gc.collect()
            ######################################################

            ##################### Cliques of size 3
            allTissueTriangles = [[A,B,C] for [A,B,C,bit] in U["chrValues"][ch]["cliques"] if ((bit&tisBit)==tisBit)]
            nodesInAllCliques = set([c[i] for i in [0,1,2] for c in allTissueTriangles])
            linksInAllCliques = set([(c[p[0]], c[p[1]]) for p in [[0,1], [0,2], [1,2]] for c in allTissueTriangles])
            header+=["Clique count", "Unique nodes in cliques", "Unique links in cliques"]
            row+=[len(allTissueTriangles), len(nodesInAllCliques), len(linksInAllCliques)]
            
            curStats = getExpressionsForSegments(nodesInAllCliques, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3segs {typ}: for all segments that are in any clique"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("Clique nodes", []).append(curStats)
            del nodesInAllCliques
            del linksInAllCliques


            #####################Cliques where all 3 nodes have an encode property 
            encS = U["chrValues"][ch]["encodeStates"]
            tissueTriangles = [c3 for c3 in allTissueTriangles if encS[c3[0]].get(tis,0)>0 and encS[c3[1]].get(tis,0)>0 and encS[c3[2]].get(tis,0)>0]
            del allTissueTriangles
            nodesInCliques = set([c[i] for i in [0,1,2] for c in tissueTriangles])
            linksInCliques = set([(c[p[0]], c[p[1]]) for p in [[0,1], [0,2], [1,2]] for c in tissueTriangles])
            header+=["CCRe clique count: cliques where all 3 nodes have a CCRe property", 
                    "Unique CCRe C3 nodes: node count in cliques where all 3 nodes have a CCRe property", 
                    "Unique CCRe C3 links: link count in cliques where all 3 nodes have a CCRe property"]
            row+=[len(tissueTriangles), len(nodesInCliques), len(linksInCliques)]
            
            curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3EncS {typ}: for all segments that are in a clique where all 3 nodes have some encode prop"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("Encode clique nodes", []).append(curStats)

            #####################How many nodes with property P in clique nodes
            cliqueNodesWithP = [node for node in nodesInCliques if node in segmentsWithP]
            cliqueNodesWithE = [node for node in nodesInCliques if node in segmentsWithE]
            cliqueNodesWithC = [node for node in nodesInCliques if node in segmentsWithC]
            cliqueNodesWithD = [node for node in nodesInCliques if node in segmentsWithD]

            header+=["P nodes in C3: promoter segments in clique nodes", 
                    "E nodes in C3: distal enhancer segments in clique nodes", 
                    "C nodes in C3: CTCF segments in clique nodes",
                    "D nodes in C3: DNase-only segments in clique nodes"]
            row+=[len(cliqueNodesWithP), len(cliqueNodesWithE), len(cliqueNodesWithC), len(cliqueNodesWithD)]
            
            

            curStats = getExpressionsForSegments(cliqueNodesWithP, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3_P {typ}: clique nodes with a promoter in a set of cliques where all 3 nodes have some encode prop"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("promoter C3* nodes", []).append(curStats)
                
            curStats = getExpressionsForSegments(cliqueNodesWithE, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3_E {typ}: clique nodes with a distal enhancer in a set of cliques where all 3 nodes have some encode prop"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("Distal enhancer C3* nodes", []).append(curStats)

            curStats = getExpressionsForSegments(cliqueNodesWithC, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3_C {typ}: clique nodes with a CTCF in a set of cliques where all 3 nodes have some encode prop"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("CTCF C3* nodes", []).append(curStats)

            curStats = getExpressionsForSegments(cliqueNodesWithD, ch, tis)
            for typ in ["GTEx", "FANTOM5"]:
                for prop in statisticsProps:
                    header+=[f"{prop} C3_D {typ}: clique nodes with a DNase-only in a set of cliques where all 3 nodes have some encode prop"]
                    row+=[curStats[typ][prop]]
                    
            R.setdefault("DNase-only C3* nodes", []).append(curStats)


            cliquesWithProperty = {}
            for triplet in combinationsPECD:
                tripletStr = getTripletStr(triplet) 
                header+=[f"{tripletStr}: clique count where the {len(triplet)} nodes these properties"]
                cliquesWithProperty[tripletStr] = [clique for clique in tissueTriangles if cliqueSatisfiesEncodeTriplet(clique, ch, tis, triplet)]
                row+=[len(cliquesWithProperty[tripletStr])]
                

                curStats = getExpressionsForListOfCliques(cliquesWithProperty[tripletStr], ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{tripletStr} {prop} {typ}"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"{tripletStr} nodes", []).append(curStats)

            

            #################### Calculate stats for cliques of size 3, where nodes have specific combinations of cCRE properties
            #E.g. PP* means 2 nodes have promoters, and the third has one, other property
            for letter in letters:
                #letters = ['P', 'E', 'C', 'D']
                EEstar = f"{letter}{letter}*"
                EE_ = [triplet for triplet in cliquesWithProperty.keys() if triplet.count(letter)>=2]
                cliquesWithProperty[EEstar] = []
                for eekey in EE_:
                    cliquesWithProperty[EEstar].extend(cliquesWithProperty[eekey])
                iii=0
                cliquesWithProperty[EEstar] = list(set([tuple(el) for el in cliquesWithProperty[EEstar]]))
                iii=0

                header+=[f"{EEstar}: clique count where the 2 nodes have these properties and the third has any other property"]
                row+=[len(cliquesWithProperty[EEstar])]
                curStats = getExpressionsForListOfCliques(cliquesWithProperty[EEstar], ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{EEstar} {prop} {typ}"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"{EEstar} nodes", []).append(curStats)

            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            addCliquesAndAlmostCliques(U, ch, tis, 5) #fills only with data for one tissue
            #Calculated cliques of size 4 and 5
            for csize in [4,5]:
                ##################################################################
                #Process cliques of size 4 --> csize and more

                #Find all cliques4 that are in the current tissue
                # allTissueCliques4 = [[A,B,C,D] for [A,B,C,D,bit] in U["chrValues"][ch]["cliques4"] if ((bit&tisBit)==tisBit)]
                # Extract cliques of the specified size (csize) and filter based on the tisBit condition

                allTissueCliques4 = [[*clique[:csize]] for clique in U["chrValues"][ch][f"cliques{csize}"] if (clique[-1] & tisBit) == tisBit]

                #////
                nodesInAllCliques_tmp = set([c[i] for i in range(csize) for c in allTissueCliques4])
                linksInAllCliques_tmp = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in allTissueCliques4])
                header+=[f"Clique{csize} count", f"Unique nodes in cliques{csize}", f"Unique links in cliques{csize}"]
                row+=[len(allTissueCliques4), len(nodesInAllCliques_tmp), len(linksInAllCliques_tmp)]
                curStats = getExpressionsForSegments(nodesInAllCliques_tmp, ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{prop} C{csize}segs {typ}: for all segments that are in any clique{csize}"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"Clique{csize} nodes", []).append(curStats)
                #////


                #Cliques4 where all nodes have some cCRE property
                # tissueCliques4 = [c4 for c4 in allTissueCliques4 if encS[c4[0]].get(tis,0)>0 and encS[c4[1]].get(tis,0)>0 and encS[c4[2]].get(tis,0)>0 and encS[c4[3]].get(tis,0)>0]
                # Filter cliques where each element's corresponding encS value for tis is greater than 0
                tissueCliques4 = [c4 for c4 in allTissueCliques4 if all(encS[node].get(tis, 0) > 0 for node in c4)]
                #maybe calculate the nodes C4 again? 

                #Collect basic statistics about C4 cliques
                # nodesInCliques = set([c[i] for i in [0,1,2,3] for c in tissueCliques4])
                nodesInCliques = {node for c in tissueCliques4 for node in c}

                linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in tissueCliques4])
                
                header+=[f"CCRe clique{csize} count: cliques{csize} where all {csize} nodes have a CCRe property", 
                        f"Unique CCRe C{csize} nodes: node count in cliques{csize} where all {csize} nodes have a CCRe property", 
                        f"Unique CCRe C{csize} links: link count in cliques{csize} where all {csize} nodes have a CCRe property"]
                row+=[len(tissueCliques4), len(nodesInCliques), len(linksInCliques)]
                
                curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{prop} C{csize}EncS {typ}: for all segments that are in a clique{csize} where all {csize} nodes have some encode prop"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"Encode clique{csize} nodes", []).append(curStats)


                #Look at CCCC, PPPP, DDDD, EEEE cliques
                #use segmentsWithD, segmentsWithC, segmentsWithP, segmentsWithE
                #Those are sets with all nodes that have this property
                #use tissueCliques4

                for letter in ["C", "P", "D", "E"]:
                    XXXX = letter*csize #E.g. CCCC or PPPPP
                    segmentsWithX=None
                    if letter=="C": segmentsWithX = segmentsWithC
                    elif letter=="P": segmentsWithX = segmentsWithP
                    elif letter=="D": segmentsWithX = segmentsWithD
                    elif letter=="E": segmentsWithX = segmentsWithE
                    
                    cliques4WithXXXX = [cli for cli in tissueCliques4  if len([cli_n for cli_n in cli if cli_n in segmentsWithX])==csize]
                    nodesInCliques = set([c[i] for i in range(csize) for c in cliques4WithXXXX])
                    linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in cliques4WithXXXX])
                    header+=[f"{XXXX} CCRe clique{csize} count: cliques{csize} where all {csize} nodes have this CCRe property", 
                            f"{XXXX} Unique CCRe C{csize} nodes: node count in cliques{csize} where all {csize} nodes have this CCRe property", 
                            f"{XXXX} Unique CCRe C{csize} links: link count in cliques{csize} where all {csize} nodes have this CCRe property"]
                    row+=[len(cliques4WithXXXX), len(nodesInCliques), len(linksInCliques)]
                    
                    curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                    for typ in ["GTEx", "FANTOM5"]:
                        for prop in statisticsProps:
                            header+=[f"{prop} {XXXX} {typ}: for all segments that are in a clique{csize} where all {csize} nodes have these encode props"]
                            row+=[curStats[typ][prop]]
                            
                    R.setdefault(f"{XXXX} nodes in C{csize}", []).append(curStats)
                

                #Calculate PP**, EE**, DD**, CC**
                def howManyObjectsWithThisProperty(obj, setOfNodesWithDesiredProperty):
                    #obj is a list, an iterable
                    return len([el for el in obj if el in setOfNodesWithDesiredProperty])
                    
                for letter in ["C", "P", "D", "E"]:
                    XXXX = letter*(csize-2)+"**"
                    segmentsWithX=None
                    if letter=="C": segmentsWithX = segmentsWithC
                    elif letter=="P": segmentsWithX = segmentsWithP
                    elif letter=="D": segmentsWithX = segmentsWithD
                    elif letter=="E": segmentsWithX = segmentsWithE
                    
                    cliques4WithXXXX = [cli for cli in tissueCliques4 if howManyObjectsWithThisProperty(cli[:-1], segmentsWithX)>=(csize-2)]
                    nodesInCliques = set([c[i] for i in range(csize) for c in cliques4WithXXXX])
                    linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in cliques4WithXXXX])
                    header+=[f"{XXXX} CCRe clique{csize} count: cliques{csize} where the nodes have these CCRe properties", 
                            f"{XXXX} Unique CCRe C{csize} nodes: node count in cliques{csize} where the nodes have these CCRe properties", 
                            f"{XXXX} Unique CCRe C{csize} links: link count in cliques{csize} where the nodes have these CCRe properties"]
                    row+=[len(cliques4WithXXXX), len(nodesInCliques), len(linksInCliques)]
                    
                    curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                    for typ in ["GTEx", "FANTOM5"]:
                        for prop in statisticsProps:
                            header+=[f"{prop} {XXXX} {typ}: for all segments that are in a clique{csize} where the {csize} nodes have these encode props"]
                            row+=[curStats[typ][prop]]
                            
                    R.setdefault(f"{XXXX} nodes in C{csize}", []).append(curStats)    
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #Process almost cliques
            postfix="_minus1"
            for csize in [4,5]:
                ##################################################################
                #Process almost cliques of size 4 --> csize and more
                allTissueCliques4 = [[*clique[:csize]] for clique in U["chrValues"][ch][f"cliques{csize}_minus1"] if (clique[-1] & tisBit) == tisBit]+[[*clique[:csize]] for clique in U["chrValues"][ch][f"cliques{csize}"] if (clique[-1] & tisBit) == tisBit]
                allTissueCliques4 = [tuple(cli) for cli in allTissueCliques4]
                allTissueCliques4 = list(set(allTissueCliques4))
                #////
                nodesInAllCliques_tmp = set([c[i] for i in range(csize) for c in allTissueCliques4])
                linksInAllCliques_tmp = set()
                for p in list(itertools.combinations(list(range(0,csize)), 2)):
                    for c in allTissueCliques4:
                        linksInAllCliques_tmp.add(tuple([c[p[0]], c[p[1]]]))
                header+=[f"Clique{csize}{postfix} count", f"Unique nodes in cliques{csize}{postfix}", f"Unique links in cliques{csize}{postfix}"]
                row+=[len(allTissueCliques4), len(nodesInAllCliques_tmp), len(linksInAllCliques_tmp)]
                
                curStats = getExpressionsForSegments(nodesInAllCliques_tmp, ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{prop} C{csize}{postfix}segs {typ}: for all segments that are in any clique{csize}{postfix}"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"Clique{csize}{postfix} nodes", []).append(curStats)
                #////
                
                
                
                tissueCliques4 = [c4 for c4 in allTissueCliques4 if all(encS[node].get(tis, 0) > 0 for node in c4)]
                nodesInCliques = {node for c in tissueCliques4 for node in c}
                linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in tissueCliques4])
                
                header+=[f"CCRe clique{csize}{postfix} count: cliques{csize}{postfix} where all {csize} nodes have a CCRe property", 
                        f"Unique CCRe C{csize}{postfix} nodes: node count in cliques{csize}{postfix} where all {csize} nodes have a CCRe property", 
                        f"Unique CCRe C{csize}{postfix} links: link count in cliques{csize}{postfix} where all {csize} nodes have a CCRe property"]
                row+=[len(tissueCliques4), len(nodesInCliques), len(linksInCliques)]
                
                curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{prop} C{csize}{postfix}EncS {typ}: for all segments that are in a clique{csize}{postfix} where all {csize} nodes have some encode prop"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"Encode clique{csize}{postfix} nodes", []).append(curStats)


                #Look at CCCC, PPPP, DDDD, EEEE cliques
                #use segmentsWithD, segmentsWithC, segmentsWithP, segmentsWithE
                #Those are sets with all nodes that have this property
                #use tissueCliques4

                for letter in ["C", "P", "D", "E"]:
                    XXXX = letter*csize
                    segmentsWithX=None
                    if letter=="C": segmentsWithX = segmentsWithC
                    elif letter=="P": segmentsWithX = segmentsWithP
                    elif letter=="D": segmentsWithX = segmentsWithD
                    elif letter=="E": segmentsWithX = segmentsWithE
                    
                    cliques4WithXXXX = [cli for cli in tissueCliques4  if len([cli_n for cli_n in cli if cli_n in segmentsWithX])==csize]
                    nodesInCliques = set([c[i] for i in range(csize) for c in cliques4WithXXXX])
                    linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in cliques4WithXXXX])
                    header+=[f"{XXXX} CCRe clique{csize}{postfix} count: cliques{csize}{postfix} where all {csize} nodes have this CCRe property", 
                            f"{XXXX} Unique CCRe C{csize}{postfix} nodes: node count in cliques{csize}{postfix} where all {csize} nodes have this CCRe property", 
                            f"{XXXX} Unique CCRe C{csize}{postfix} links: link count in cliques{csize}{postfix} where all {csize} nodes have this CCRe property"]
                    row+=[len(cliques4WithXXXX), len(nodesInCliques), len(linksInCliques)]
                    
                    curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                    for typ in ["GTEx", "FANTOM5"]:
                        for prop in statisticsProps:
                            header+=[f"{prop} {XXXX} {typ}: for all segments that are in a clique{csize}{postfix} where all {csize} nodes have these encode props"]
                            row+=[curStats[typ][prop]]
                            
                    R.setdefault(f"{XXXX} nodes in C{csize}{postfix}", []).append(curStats)
                

                #Calculate PP**, EE**, DD**, CC**
                def howManyObjectsWithThisProperty(obj, setOfNodesWithDesiredProperty):
                    #obj is a list, an iterable
                    return len([el for el in obj if el in setOfNodesWithDesiredProperty])
                    
                for letter in ["C", "P", "D", "E"]:
                    XXXX = letter*(csize-2)+"**"
                    segmentsWithX=None
                    if letter=="C": segmentsWithX = segmentsWithC
                    elif letter=="P": segmentsWithX = segmentsWithP
                    elif letter=="D": segmentsWithX = segmentsWithD
                    elif letter=="E": segmentsWithX = segmentsWithE
                    
                    cliques4WithXXXX = [cli for cli in tissueCliques4 if howManyObjectsWithThisProperty(cli[:-1], segmentsWithX)>=(csize-2)]
                    nodesInCliques = set([c[i] for i in range(csize) for c in cliques4WithXXXX])
                    linksInCliques = set([(c[p[0]], c[p[1]]) for p in list(itertools.combinations(list(range(0,csize)), 2)) for c in cliques4WithXXXX])
                    header+=[f"{XXXX} CCRe clique{csize}{postfix} count: cliques{csize}{postfix} where the nodes have these CCRe properties", 
                            f"{XXXX} Unique CCRe C{csize}{postfix} nodes: node count in cliques{csize}{postfix} where the nodes have these CCRe properties", 
                            f"{XXXX} Unique CCRe C{csize}{postfix} links: link count in cliques{csize}{postfix} where the nodes have these CCRe properties"]
                    row+=[len(cliques4WithXXXX), len(nodesInCliques), len(linksInCliques)]
                    
                    curStats = getExpressionsForSegments(nodesInCliques, ch, tis)
                    for typ in ["GTEx", "FANTOM5"]:
                        for prop in statisticsProps:
                            header+=[f"{prop} {XXXX} {typ}: for all segments that are in a clique{csize}{postfix} where the {csize} nodes have these encode props"]
                            row+=[curStats[typ][prop]]
                            
                    R.setdefault(f"{XXXX} nodes in C{csize}{postfix}", []).append(curStats)    
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            #use tissueLinks
            for letter in ["C", "P", "D", "E"]:
                XX = letter*2+"_links"
                segmentsWithX=None
                if letter=="C": segmentsWithX = segmentsWithC
                elif letter=="P": segmentsWithX = segmentsWithP
                elif letter=="D": segmentsWithX = segmentsWithD
                elif letter=="E": segmentsWithX = segmentsWithE

                activeLinks = [link for link in tissueLinks if howManyObjectsWithThisProperty(link[:2], segmentsWithX)==2]
                nodesInLinks = set([c[i] for i in [0,1] for c in activeLinks])
                header+=[f"{XX} link count where both nodes have these cCRE property", 
                        f"{XX} node count in links where both nodes have these cCRE property"]
                row+=[len(activeLinks), len(nodesInLinks)]
                
                curStats = getExpressionsForSegments(nodesInLinks, ch, tis)
                for typ in ["GTEx", "FANTOM5"]:
                    for prop in statisticsProps:
                        header+=[f"{prop} {XX} {typ}: for all segments that are in any link where both nodes have these encode props"]
                        row+=[curStats[typ][prop]]
                        
                R.setdefault(f"{XX} nodes", []).append(curStats)  

            iii=0


            L.append(row)
    U["chrValues"][ch].clear()
    return [header, L]




def appendToFile(LL, header, resFn):
    # Check if the file exists
    file_exists = os.path.isfile(resFn)
    
    # Open the file in append mode, create it if it doesn't exist
    with open(resFn, 'a', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        # Optionally, if the file doesn't exist, write a header
        if not file_exists:
            # Writing a header
            #But before that change the header elements slightly, using 
            if os.path.exists(headerMappingFN):
                # Read the mapping from the JSON file
               
                with open(headerMappingFN, 'r') as json_file:
                    header_mapping = json.load(json_file)
                header = [header_mapping.get(headerElem, headerElem) for headerElem in header]

            spamwriter.writerow(header)
        # Write the provided list of lists to the file
        spamwriter.writerows(LL)



for chh in U["chrNames"]:
    [header, L] = calculatecCRE_stats([chh])
    appendToFile(L, header, rezFn)
    gc.collect()

