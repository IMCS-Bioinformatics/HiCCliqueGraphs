import sys
sys.stdout = open('outputRandomizationLogs.log', 'w')
sys.stderr = sys.stdout

import os
from pathlib import Path
import random
import math
import time
import copy
import json
import Utils
from Utils import LengthGroups
from Utils import SegmentGroups 
# -*- coding: utf-8 -*-
"""
UniversalDS + ChrData classes
universal.py


@author: Andr
"""
import json
import sys
import csv
import networkx as nx
import matplotlib.pyplot as plt
import sys


def getBitCount(bits): #returns number of 1s in binary form
    count = 0
    while (bits):
        bits &= (bits-1)
        count+= 1
    return count



class UniversalDS:
    # Class that holds data from the Universal file and has Hi-C data
    # It has data for all chromosomes and for all tissue types
    def __init__(self, fn, assertion=False):
        print("Initializing U")
        self.fn = fn
        print(fn)
        self.readData() #sets self.data, self.{chrs, tissueBits, bitToTisname}
        self.assertion = assertion #False by default. If true, validates that received data is of correct format
    
    def readData(self):
        print("readData called")
        f = open(self.fn)
        print("opened", self.fn)
        self.data = json.load(f) #list of links and segments
        f.close()
        self.chrs = self.data["chrNames"]
        self.tissueBits = self.data["tissueBits"] #{"aCD4": 1, "EP": 2, "Ery": 4, "FoeT": 8,...
        self.bitToTisname = { self.tissueBits[tisName]: tisName for tisName in self.tissueBits.keys() }
        self.tissues = list(self.tissueBits.keys())
        
        #set self.DS
        if self.data["pvalue"]==5: 
            self.DS="BloodCellPCHiC"
        elif self.data["pvalue"] in [6,10]: 
            self.DS="NormalHi-C"
        elif self.data["pvalue"]==0.7: 
            self.DS="pcHi-C"
        else: self.DS = "Unknown"
        #if self.fn[-9:]=="rand.json": self.DS+="-rnd"
        if "rand" in self.fn: self.DS+="-rnd"
        
    
    def getTissueNameList(self, bits):
        #returns list of tissues from a bitmap of tissue types
        return [self.data["tissueIDs"][i] 
            for i in range(len(self.data["tissueIDs"])) if (((1<<i)&bits)>0) ]
    
    def getTisPairs(self, includeSelfPairs=True):
        #makes a list of tissue pairs, includeSelfPairs:=True if (A,A) is a valid pair
        # #returns [["tis1", "tis2"], [], [], ...] 
        L = []
        t = self.tissues

        for tis1 in t:
            for tis2 in t:
                if t.index(tis2)<t.index(tis1): continue
                if not includeSelfPairs:
                    if tis1==tis2: continue
                L.append([tis1, tis2])
        return L
    
    def saveData(self, fn):
        with open(fn, 'w') as json_file:
            json.dump(self.data, json_file)
    
    
        

class ChrData:
    # Holds data of Hi-C graph for one chromosome
    def __init__(self, owner, ch, allLinks=None, minLinkTissueCount=1, tissueMask=0):
        self.owner = owner # UniverslaDS instance
        self.DS = owner.DS
        self.ch = ch # chromosome from UniversalDS
        self.segments = self.owner.data["chrValues"][ch]["segments"]
        self.minLinkTissueCount = minLinkTissueCount # criteria for link filtering - at least this many tissues must be in each link
        self.tissueMask = tissueMask# criteria for link filtering - these tissues must be in each link. Extra tissues ar aceptable.
        if type(tissueMask)==str: 
            try:
                self.tissueMask = self.owner.tissueBits[tissueMask]
                #maybe tissue name was given
            except KeyError:
                self.error("tissueMask is invalid", tissueMask)
        elif type(tissueMask)==list:
            try:
                m=0
                for tisName in tissueMask:
                    m = (m | (self.owner.tissueBits[tisName]))
                self.tissueMask = m
                #maybe list of tissue names was given
            except KeyError:
                self.error("tissueMask is invalid", tissueMask)
        #allLinks - list of links. Can be passed by user or links from UniversalDS will be used
        if allLinks is None:
            self.allLinks = self.owner.data["chrValues"][self.ch]["links"]
        else:
            self.allLinks = allLinks
        self.links = self.filterLinks() #sets self.links, based on criteria in self.minLinkTissueCount and self.tissueMask
        self.updateFunctions = []
        self.segmentIndToMidpoint = {ind: (self.segments[ind][0]+self.segments[ind][1])//2 for ind in range(len(self.segments))} #dict to get segment midpoint from segment index
    def __copy__(self):
        return ChrData(self.owner, self.ch, self.links, self.minLinkTissueCount, self.tissueMask)
    
    def warning(self, text, values):
        print("\n#### WARNING: ", text, values)
    def error(self, text, values):
        print("\n#### ERROR: ", text, values)
        sys.exit()
    
    def filterLinks(self):
        # create list of links that are in at least minTissueCount tissues
        # this method removes the dict with pvalues for each link. Removed because they are never used
        # result is a list of [Aind, Bind, bits] for current chr
        linksAll = self.allLinks
        if self.minLinkTissueCount<2 and self.tissueMask==0: return [link[:3] for link in linksAll] 
        return [link[:3] for link in linksAll if ( (getBitCount(link[2])>=self.minLinkTissueCount) and ( (link[2]&self.tissueMask)==self.tissueMask) )]
    
    def updateLinks(self, newCondition): 
        #Receives a function and keeps only those links that pass the test. It gets a link in form [A, B, bit]
        # it modifies self.links by performin filters on it
        #check if function is of valid form
        try:
            res = newCondition([10,11,2])
            if type(res)!=type(True): self.error("ChrData.updateLinks got a function that returns non-bool value")
        except TypeError:
            self.error("ChrData.updateLinks got a function that does not accept links in their current form. UpdateLinks not executed", "")
        self.links = [link for link in self.links if newCondition(link)]
        self.updateFunctions.append(newCondition)
    
    def makeAdjAndBitsOfLink(self):
        # creates an adj structure for all links and creates a dictionary bitsOfLink
        # bitsOfLink[(A,B)] == the bitmap of the link
        adj = dict()
        bitsOfLink = dict()
        for link in self.links:
            [A, B, bitmap] = link
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
    
    def getC3(self):
        # finds all triangles in current chr. At least one tissue must be shared among each link of the triangle
        # result is a list of tuples (A, B, C, bit), where A, B, C are sorted asc. and where bit is the max bitmap shared among all links (AB)&(AC)&(BC)
        # only those triangles are kept that have at least one common tissue in all 3 links
        # uses links == self.links
        # links argumentthat is used here has list of links in form [A, B, bitmap] for one chromosome
        
        [adj, bitsOfLinks] = self.makeAdjAndBitsOfLink() 
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
    
    def getAllSegments(self):
        #gets list of all segments that are a part of at least one link
        s = set()
        for li in self.links:
            for i in [0,1]:
                s.add(li[i])
        return sorted(list(s))
    
    def listOfIndsToListOfLoci(self, L):
        #uses self.segmentIndToMidpoint to do the translation
        try:
            LL = [self.segmentIndToMidpoint[el] for el in L]
        except KeyError:
            self.owner.error("KeyError in listOfIndsToListOfLoci", L)
        
        return LL
    
    def getAllSegmentsLoci(self):
        #returns list of all segments in form of loci
        segInds = self.getAllSegments()
        return self.listOfIndsToListOfLoci(segInds)
    
    def getListOfSegmentLociSegments(self):
        #returns list of nodes, each node in form [A, B] where A and B are loci of the segment
        segInds = self.getAllSegments()
        return [self.segments[i] for i in segInds]
    
    def summarize(self):
        #method summarizes data about this chromosome and returns a result in form of a dict
        rez = dict()
        #number of links and nodes total for default pvalue

        rez["Chromosome"] = self.ch
        rez["Nodes count"] = len(self.getAllSegments())
        rez["Link count"] = len(self.links)
        
        #size of the giant component
        links = [(el[0], el[1]) for el in self.links ]
        G = nx.Graph()
        G.add_edges_from(links)
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        largestCC = G.subgraph(Gcc[0])
        rez["Giant component nodes"] = len(largestCC.nodes)
        rez["Giant component links"] = len(largestCC.edges)
        rez["Link percentage in Giant component"] = ((rez["Giant component links"])/(rez["Link count"]))
        
        return rez
    
    def getLinksByTissue(self, allLinks=None):
        #for all tissues, finds the number of links
        #allLinks - these links are used. Default - all filtered links.
        #returns {tis1: 156, tis2: 245, ...} with link counts for each tissue
        if allLinks is None: allLinks = self.links

        U = self.owner
        counts = {tis: 0 for tis in U.tissues}

        for link in allLinks:
            [A,B,bit] = link
            curTissues = U.getTissueNameList(bit)
            for tis in curTissues:
                counts[tis]+=1
        
        return counts
    
    def getNodeDegrees(self):
        #returns dict with keys in segment IDs, values - ints. {seg1: deg1, seg2: deg2}
        degs = {}
        for link in self.links:
            [A,B,bit] = link
            if A not in degs: degs[A]=1
            else: degs[A]+=1

            if B not in degs: degs[B]=1
            else: degs[B]+=1

        return degs
    
    def getLinkLength(self, link):
        #link is in form [A,B,bits,{}], returns segmentMid[B]-segmentMid[A]
        return (self.segmentIndToMidpoint[link[1]] - self.segmentIndToMidpoint[link[0]])
        
        
        
        
        
    


def getNthElementFromSet(S, n):
    #S - set; n- element index
    #This iterates over a set and gets the n-th value by hash. It an be considered a random element from set if n is chosen randomly
    i=0
    if n>=len(S):
        raise Exception("Can not get {}th element from a {}-element set".format(n, len(S)))
    for el in S:
        if i==n: break
        i+=1
    return el

class NewSmartRandomizer:
    def __init__(self, U, randomizablePart=0.80, resultFn="rezz/allBloodRandDefault.json", chrs = None, assertions=False):
        self.U = U
        if chrs is None: chrs=list(reversed(U.chrs))
        self.randomizablePart=randomizablePart
        #random.seed(12345)
        #random.seed()
        self.assertions = assertions
        self.maxLenEpsilon = 0.1 #epsilon for allowed links over 2e6

        self.countersToFile = dict()

        # for ch in chrs:
        #chrs = ['chr12']
        for ch in chrs:
            self.randomizeChr(ch)
        
            #exit()

            self.saveFile(resultFn)
            self.chrCounters["final"] = copy.deepcopy(self.counters)
            self.countersToFile[ch] = copy.deepcopy(self.chrCounters)
            
            json_object = json.dumps(self.countersToFile, indent=4)
            # Writing to sample.json
            with open("{}-counters.json".format(resultFn[:-5]), "w") as outfile:
                outfile.write(json_object)

    
    def randomizeChr(self, ch):
        self.initializeChr(ch)
        init_start_time = self.start_time
        self.startSwaps()
        self.links+=self.linksTooLong
        self.links = sorted(self.links, key=lambda x: (x[0],x[1]))

        self.printCounters("After chromosome is processed")
        self.chrCounters["afterSwaps"] = copy.deepcopy(self.counters) # kas tas ir????
        self.start_time = init_start_time
        print('\n===> Total swap time:', self.getTime(), '\n')

    def getTime(self):
        ctime = time.time() - self.start_time
        self.start_time = time.time()
        return ctime
    
    def initializeChr(self, ch):
        chData = ChrData(self.U, ch)
        self.ch = ch
        self.chData = chData
        self.links = chData.allLinks #list of [A,B,bit,{}]
        self.linksToRandomize = int(len(self.links)* self.randomizablePart) #how many links to randmize
        self.linksWithWrongLength=0
        #Ņem pusi no nerandomizējamiem linkiem, un vienkārši izmet tos ārā
        allLinkLengths = []
        for link in self.links:
            allLinkLengths.append(self.chData.getLinkLength(link))
        allLinkLengths = sorted(allLinkLengths, reverse=True)
        ind = round(0.025*len(allLinkLengths))
        maxAllowedlength = allLinkLengths[ind]
        linksToKeep = []
        self.linksTooLong = []
        for link in self.links:
            curLen = self.chData.getLinkLength(link)
            if curLen > maxAllowedlength:
                self.linksTooLong.append(link)
            else:
                linksToKeep.append(link)
        
        self.links = linksToKeep


        self.maxAllowedLength = max([self.getLinkLength(link) for link in self.links])

        self.linkCount = len(self.links)
        self.segmentCount = len(self.chData.segments)
        self.chrCounters = {}

        self.lengthGroups = LengthGroups(self, self.links, minElementCount=32, minLengthCount=2)

        # self.delta = self.linkCount >> 11 # == // 1024
        self.delta = int(math.log2(self.linkCount))
        # self.neighbourGroups = None
        self.start_time = time.time()

        
        # self.shakeLinksWithWrongLength=0

        print("Randomizing {ch}\n\ntotal links: {tl}". format(ch=ch, tl=len(self.chData.links)))
        print("{} will be randomized".format(self.linksToRandomize))
        print("Segment count: {};\n".format(len(self.chData.segments)))

        self.changedLinkInds = set() #set of link indeces

        # self.unchangedLinkInds = [i for i in range(self.linkCount)]
        # random.shuffle(self.unchangedLinkInds)
        # #self.unchangedLinkIndMap = {i: self.unchangedLinkInds[i] for i in range(self.linkCount)}
        # self.unchangedLinkIndMap = {self.unchangedLinkInds[i] : i for i in range(self.linkCount)}
        # self.unchangedLinkInds = [i for i in range(self.linkCount)]
        unchangedLinkInds = [i for i in range(self.linkCount)]
        random.shuffle(unchangedLinkInds)
        self.unchangedLinks = Utils.DinamicList(unchangedLinkInds, 'unchangedLinks')


        self.linksWithNoCandidates = set()

        self.linkIndInList = [i for i in range(len(self.links))] # link ind in unchangedLinkInds...
        self.linkIndType = [0 for i in range(len(self.links))] # 0: unchanged, 1: changed, 2: no candidates
        self.originalLinks = set([(el[0],el[1]) for el in self.links]) #set of (A,B)
        self.existingLinks = set([(el[0],el[1]) for el in self.links]) #set of (A,B)

        allSegments = set([link[i] for link in self.links for i in [0,1]])
        maxSeg = max(allSegments)
        self.linksAt = {}
        segmConnCounts = [0 for _  in range(maxSeg+1)]
        for link in self.links:
            segmConnCounts[link[0]]+=1
            segmConnCounts[link[1]]+=1
        iii=9
        minSegmentCountParam = round(math.log2(len(allSegments)))//4 #+-32
        minElementCountParam = round(math.log2(self.linkCount)*minSegmentCountParam) #512
        
        print("32-", minSegmentCountParam, "and 512-", minElementCountParam)
        self.segmGroups = SegmentGroups(self.links, minElementCount=minElementCountParam, minSegmentCount=minSegmentCountParam, txt='segmGroups') #ar param var speleties

        # linkIndex = random.randrange(0, len(self.links))
        # self.getSelectedCandidateIndexNew(linkIndex)

        self.counters = {}
        self.counters = {
            "ch": self.ch,
            "time": 0,
            "totalLinks": len(self.chData.links),
            "totalSegments": len(self.chData.segments),
            "linksToRandomize": self.linksToRandomize,
        }
    def getLinkLength(self, link):
        #link is in form [A,B,bits,{}], returns segmentMid[B]-segmentMid[A]
        return self.linkLengthFromPoints(link[0], link[1])
    
    def linkLength(self, linkInd):
        return self.linkLengthFromPoints(self.links[linkInd][1], self.links[linkInd][0])
    
    def linkLengthFromPoints(self, A, B):
        #input - 2 points (2 indeces of segments (link endpoints))
        #output  -AB length on the chromosome, in BP
        return abs(self.chData.segmentIndToMidpoint[B] - self.chData.segmentIndToMidpoint[A])
    
    def normLinkLength(self,  linkInd):
        #paņem linku, paņem abu galu reālās koordinātas, dalīs ar 2048
        return self.normLength(self.linkLength(linkInd))
    def normLength(self, length): #paņem garumu, normē to
        return length>>12
    
    def getLinksWithWrongSize(self):
        return sum([el for el in self.boxWeights.values() if el>0])*2
        
    def printCounters(self, text="Unknown context"):
        self.counters["time"] = self.getTime()
        self.counters["SwappedLinkCount"] = len(self.changedLinkInds)
        self.counters["Unchanged links"] = self.unchangedLinks.length()
        self.counters["LinksWithWrongSize"] = self.lengthGroups.badLinkLengthCount
        #self.groupInitLength = {g: len(self.groupElem[g]) for g in self.groupElem.keys()}
        c = 0
        for gr in self.lengthGroups.groupElem.keys():
            c+=abs(len(self.lengthGroups.groupElem[gr])-self.lengthGroups.groupInitLength[gr])
        self.counters["calculated Wrong sized links"] = c
        self.counters["text"] = text

        with open("RandomizationLogs.txt", "a") as file_object:
            file_object.write(str(self.counters))
            file_object.write("\n")


        print(text)
        print(self.counters)
        print("\n\n")
    def startSwaps(self):
        #galvenais
        prevPrint=0
        self.printCounters("First print")
        #swapStep = self.linkCount // self.delta #delta==log(links)
        while len(self.changedLinkInds) < self.linksToRandomize and self.unchangedLinks.length()>0 :
            swaps_start_time = self.start_time
            swapStep = int((0.05 + random.random()*0.05)*self.linksToRandomize)
            # swapStep = 1000
            self.swapUnchecked(swapStep) #gribam tik daudz apmest
            startSwaps_start_time = swaps_start_time
            print('\n===> swapUnchecked is finished')
            self.printCounters("After throw")
            self.chrCounters["After throw"] = copy.deepcopy(self.counters)
            # self.assertMapAndListValidity()
            self.postProcess()
            #self.start_time = swaps_start_time
            print('\n===> spostProcess finished:')
            self.printCounters("After throw + shake")
            self.chrCounters["After throw + shake"] = copy.deepcopy(self.counters)
            for _ in range(3):
                print("---------------------------------------------------------------")
            iii = 1

    def getCandidateScore(self, testLinkInd, candLinkInd):
        # (U,V,A,B)
        [benefitScore1, disp1, value1] = self.candidateScore([candLinkInd[0], candLinkInd[1], testLinkInd[0], testLinkInd[1]])
        # # (V,U,A,B)
        # [benefitScore2, disp2, value2] = self.candidateScore([candLinkInd[1], candLinkInd[0], testLinkInd[0], testLinkInd[1]])
        # benefitScore = max(benefitScore1, benefitScore2)#best we can do with these links
        return benefitScore1
        
    def isUndefScore(self, score):
        return score[0] == -100

    def undefScore(self):
        return [-100, sys.maxsize, [], []]

    def isUnchanged(self, linkInd):
        return self.linkIndType[linkInd] == 0

    def isChanged(self, linkInd):
        return self.linkIndType[linkInd] == 1

    def isNoCandidates(self, linkInd):
        return self.linkIndType[linkInd] == 2

    def candidateScore(self, p): # p==[A,B,V,U] creates AU and BV. returns [benefitScore, disp1,...]
        if len(set(p[:4]))<4:
            return self.undefScore()
        # [4, 242, [2082, 2124, 2071, 2029, 4951, 4885], [42, 42, 53, 53]]
        if p[4] == 4951 and p[5] == 4885:
            iii = 1

        newPoints = [tuple(sorted([p[0], p[3]])), tuple(sorted([p[1], p[2]]))]
        if newPoints[0] in self.existingLinks or newPoints[1] in self.existingLinks or \
            newPoints[0] in self.originalLinks or newPoints[1] in self.originalLinks:
            #if newPoints[0] in self.originalLinks or newPoints[1] in self.originalLinks:
            return self.undefScore()

        #lengths = [abs(p[0] - p[1]), abs(p[2] - p[3]), abs(p[1] - p[2]), abs(p[0] - p[3])]
        lengths = [self.linkLengthFromPoints(p[0],p[1]), self.linkLengthFromPoints(p[2],p[3]),self.linkLengthFromPoints(p[1],p[2]),self.linkLengthFromPoints(p[0],p[3])]
        normLengths = [self.normLength(el) for el in lengths] #normalized
        minLengths = sorted([min(lengths[:2]), min(lengths[2:])])
        #lengthGroups = [self.lengthGroups.lenGroup[v] for v in normLengths]
        lengthGroups = []
        for v in normLengths:
            if v >= len(self.lengthGroups.lenGroup):
                ###link is too long, longer than double anything there was in original graph!!
                return self.undefScore() 
            lengthGroups.append(self.lengthGroups.lenGroup[v])
            
        copyBoxWeight = {}
        for groupG in lengthGroups:
            copyBoxWeight[groupG] = self.lengthGroups.linkGroupLengthDiff(groupG) #tekosais - sakotnejais
            
        newBenefitScore=0
        for i in [-1, 1]:
            for j in [0, 1]:
                if copyBoxWeight[lengthGroups[i + j + 1]] * i < 0:
                    newBenefitScore+=1
                else:
                    newBenefitScore -= 1
                copyBoxWeight[lengthGroups[i + j + 1]] += i

        oLength = sorted([lengths[0], lengths[1]]) + sorted([lengths[2], lengths[3]])
        minRel = minLengths[0] / minLengths[1]
        if (newBenefitScore < 4 and minRel<0.25) or (oLength[3]>self.maxAllowedLength):
            #print("Eliminated long link.", lengths, self.maxAllowedLength)
            return self.undefScore()
        
        disp1 = round(((oLength[0] - oLength[2])**2 + (oLength[1]- oLength[3])**2) / \
            ((oLength[0] - oLength[1])**2 + (oLength[2]- oLength[3])**2 + 1), 4)
        # if 0 in [oLength[0], oLength[1]]:
        #     iii=9
        value = ((oLength[2] - oLength[0]) ** 2) / oLength[0] + ((oLength[3] - oLength[1]) ** 2) / oLength[1]
        disqualified = round(math.log10(1 + math.log10(1 + value)), 4)

        if newBenefitScore > 2:
            iii = 1
        newBenefitScore -= disqualified
        if disqualified < 0.3:
            iii = 1

        return [newBenefitScore, disp1, p, lengths]

    def getCoord(self, a, b):
        return [self.links[a][0], self.links[a][1], self.links[b][0], self.links[b][1]]
    

    def changeRandUnchangedLink(self, uncheckedLinkInd):
        candidates = []
        candidateLinkInd = self.lengthGroups.anyCandidate(uncheckedLinkInd)
        if candidateLinkInd<0:
            return []
        coord1 = self.getCoord(uncheckedLinkInd, candidateLinkInd) + [uncheckedLinkInd, candidateLinkInd]
        coord2 = [coord1[0], coord1[1], coord1[3], coord1[2], uncheckedLinkInd, candidateLinkInd]
        c1 = self.candidateScore(coord1)
        if not self.isUndefScore(c1):
            candidates.append(c1)
        c1 = self.candidateScore(coord2)
        if not self.isUndefScore(c1):
            candidates.append(c1)
        return candidates
        #return [self.candidateScore(coord1), self.candidateScore(coord2)]

    def assertMapAndListValidity(self):
        if not self.unchangedLinks.isValid():
            iii=9
            exit()


    def swapUnchecked(self, changedCount): #param - skaits, cik grib apmest, this is unsupervised swapping
        n = changedCount
        while n > 0 and self.unchangedLinks.length() > 0 and len(self.changedLinkInds) < self.linksToRandomize:
            #self.assertMapAndListValidity() #This fails
            uncheckedLinkInd = self.unchangedLinks.randValue()
            #self.unchangedLinks.removeValue(uncheckedLinkInd)
            # self.startNeighbours(uncheckedLinkInd, True, 1) #linka apkartnē savāc linkus, 4 ir multiplikators
            # if len(self.neighbourCandidates) > 0:
            #     #score = max(self.changeRandUnchangedLink(uncheckedLinkInd))
            #     score = sorted(self.changeRandUnchangedLink(uncheckedLinkInd), key=lambda x: (-x[0], x[1]))
            #     if not self.isUndefScore(score[0]):
            #         n -= 1
            #         self.doTheSwap(score[0])
            # n -= 1
            # score = sorted(self.changeRandUnchangedLink(uncheckedLinkInd), key=lambda x: (-x[0], x[1]))
            # chosenLinkIndex = self.changeRandUnchangedLink(uncheckedLinkInd)
            score = self.getSelectedCandidateIndexNew(uncheckedLinkInd)
            if not self.isUndefScore(score):
                self.doTheSwap(score)
            n-=1
            iii=9
        # self.assertMapAndListValidity()

    def getAllPositiveChangedLinks(self):
        L = []
        for changedLinkInd in self.changedLinkInds:
            if self.boxWeights[self.linkLength(changedLinkInd)]>0:
                L.append(changedLinkInd)
        random.shuffle(L)
        return L

    def postProcess(self):
        print("-------------------------------------------------------------------------------------")
        linksToProcess = (len(self.changedLinkInds) / self.linksToRandomize)*self.linkCount
        count=0
        T=512
        improvementPerT = max(1,round(math.log2(self.lengthGroups.badLinkLengthCount)/2))
        self.counters["pp"] = ["links to process:", linksToProcess, "improvementPerT:", improvementPerT]
        ar = [0 for _ in range(T)]
        improvedLinkCount=0
        cancelledChangeLinkCount=-100 #?
        #1) Ja pozitīvo bokšu vairāk nav, beidzam
        #2) Nevar būt visi pozitīvie bokši. Pozitīvo ir ne vairāk kā 1/2 no visiem bokšiem. Cik ir pozitīvo linku - tas ir jāzina (self.lengthGroups.badLinkLengthCount)
        #3) improvedLinkCount nevar pārsniegt badLinkCount sākumā

        
        while self.unchangedLinks.length() > 0 and len(self.lengthGroups.positiveGroups.data)>0 and (count<linksToProcess or (improvedLinkCount/(cancelledChangeLinkCount+1))>0.5) :
            # self.assertMapAndListValidity()
            ar[count%T]=self.lengthGroups.badLinkLengthCount
            count+=1
            
            if (count%T)==0:
                s=0
                for i in range(T-1):
                    s+=(ar[i]-ar[i+1])**2
                if s<improvementPerT:
                    print("h", s, improvedLinkCount)
                    break
            
            if (count%1000)==0: 
                self.counters["cancelledChangeLinkCount"] = cancelledChangeLinkCount
                #self.printCounters("postProcess after {}th iteration".format(count))
                self.chrCounters[str(count)+"duringPostPr"] = copy.deepcopy(self.counters)

                self.counters["cancelledChangeLinkCount"] = 0
                self.counters["count"] = count
                #cancelledChangeLinkCount = 0
                


            chosenLinkIndex = self.lengthGroups.positiveCandidate(-1)
            score = self.getSelectedCandidateIndexNew(chosenLinkIndex)
            if self.isUndefScore(score) or score[0]<=0: 
                cancelledChangeLinkCount+=1
                if cancelledChangeLinkCount/count>0.95:
                    print ("cancelledChangeLinkCount 95%", cancelledChangeLinkCount, count)
                    break
                continue
            self.doTheSwap(score)
            improvedLinkCount+=1
        #self.assertMapAndListValidity()
        print(count, linksToProcess, improvedLinkCount, cancelledChangeLinkCount)

        ##emptify 
        if random.random()<0.25:
            noCandidates = self.linksWithNoCandidates.difference(self.changedLinkInds)
            for linkWithoutCandidate in noCandidates:
                self.unchangedLinks.addValue(linkWithoutCandidate)
            self.linksWithNoCandidates.clear()
            print("Emptified!")

    def getSelectedCandidateIndexNew(self, chosenLinkIndex):
        around = 1
        pool = []
        isSwapped = chosenLinkIndex in self.changedLinkInds
        for i in [0, 1]:
            index = chosenLinkIndex * 2 + i
            baseGroup = self.segmGroups[index]
            for group in range(max(0,baseGroup-around), min(baseGroup+around, self.segmGroups.groupElemMaxKey)):
                if group not in self.segmGroups.groupElem:
                    print('\n###getSelectedCandidateIndexNew:: has not group', chosenLinkIndex, index, baseGroup, group, '\n')
                    self.segmGroups.groupElem[group] = []
                    iii = 1
                linkInds = copy.copy(self.segmGroups.groupElem[group])                 
                pool += linkInds

        poolSet = set(pool)

        poolSet.remove(chosenLinkIndex * 2 + 0)
        poolSet.remove(chosenLinkIndex * 2 + 1)
        pool=list(poolSet)
        if len(poolSet)<2:
            print("poolSet is empty", chosenLinkIndex, self.segmGroups[chosenLinkIndex * 2 + 0], self.segmGroups[chosenLinkIndex * 2 + 1])
        random.shuffle(pool)
        candidates = []
        N = len(pool)//2
        count = N
        while count > 0 and len(candidates)<256:
            count -= 1
            candidateLinkIndex = pool[count] >> 1
            coord1 = self.getCoord(chosenLinkIndex, candidateLinkIndex) + [chosenLinkIndex, candidateLinkIndex]
            score = self.candidateScore(coord1)
            if not self.isUndefScore(score):
                candidates.append(score)
            
            coord2 = [coord1[0], coord1[1], coord1[3], coord1[2], chosenLinkIndex, candidateLinkIndex]
            score = self.candidateScore(coord2)
            if not self.isUndefScore(score):
                candidates.append(score)
            # self.assertMapAndListValidity()
        if len(candidates)==0:
            #print("No candidates found", chosenLinkIndex)
            self.linksWithNoCandidates.add(chosenLinkIndex)
            if self.unchangedLinks.hasValue(chosenLinkIndex):
                self.unchangedLinks.removeValue(chosenLinkIndex)
            # else:
            #     print("should not happen")
            #self.unchangedLinks.removeValue(chosenLinkIndex)
        candidates = sorted(candidates, key=lambda x: (-x[0], x[1]))
        iii=9
        cN = len(candidates)
        if cN==0: return self.undefScore() #interesanti, kurai grupai pieder. Uzkrat info par to
        N = cN
        candidateLabums = {-4:0, -2:0, 0:0, 2:0, 4:0}
        for i in range(len(candidates)):
            cand = candidates[i]
            if cand[0]>2: candidateLabums[2]=i
            elif cand[0]>0: candidateLabums[0]=i
            elif cand[0]>-2: candidateLabums[-2]=i
            else: candidateLabums[-4]=i
        for i in [2,0,-2,-4]:
            if candidateLabums[i]>=8:
                N = candidateLabums[i]
                break

        selectedCandidateIndex = math.floor(((random.random())**4)*(N)) #Seit vajag, lai nostradatu randoms, jabut pietiekosi daudz
        #self.neighbourCandidates.clear()
        iii=9
        return candidates[selectedCandidateIndex]

    def getSelectedCandidateIndex(self, chosenLinkIndex):
        #First finds several candidates for a given link
        # self.startNeighbours(chosenLinkIndex, False, 4)
        #this is list of neigbours links choose candidates from this list
        # N = len(self.neighbourCandidates)
        N = self.lengthGroups.badLinkLengthCount // 2
        candidates = []
        count = N
        # while count > 0 and len(self.neighbourCandidates)>0 and len(candidates)<100:
        while count > 0 and self.lengthGroups.badLinkLengthCount>0 and len(candidates)<100:
            count -= 1
            # candidateLinkIndex = self.nextCandidate()
            candidateLinkIndex = self.lengthGroups.positiveCandidate(chosenLinkIndex)
            if (candidateLinkIndex==64509 or candidateLinkIndex==20936):
                iii=9
            # self.assertMapAndListValidity()
            coord1 = self.getCoord(chosenLinkIndex, candidateLinkIndex) + [chosenLinkIndex, candidateLinkIndex]
            score = self.candidateScore(coord1)
            if not self.isUndefScore(score):
                candidates.append(score)
            
            coord2 = [coord1[0], coord1[1], coord1[3], coord1[2], chosenLinkIndex, candidateLinkIndex]
            score = self.candidateScore(coord2)
            if not self.isUndefScore(score):
                candidates.append(score)
            # self.assertMapAndListValidity()
        
        candidates = sorted(candidates, key=lambda x: (-x[0], x[1]))
        iii=9
        N = len(candidates)
        if N==0: return self.undefScore()
        selectedCandidateIndex = math.floor(((random.random())**4)*(N//2))
        #self.neighbourCandidates.clear()
        iii=9
        return candidates[selectedCandidateIndex]

    def remove(self, index):
        self.lengthGroups.remove(index)
        self.segmGroups.remove(index)
    def add(self, index):
        self.lengthGroups.add(index)
        self.segmGroups.add(index)
        
    def doTheSwap(self, score):
        # print('doTheSwap', score)
        #score==[benefit, dispersija, [A,B,U,V,ABIndex,UVIndex], [len(AB),len(UV),len(UB), len(AV)]]
        #self.linksWithWrongLength -= score[0]
        links = score[2][4:6]  #link indeces
        if links[0]==64509 and links[1]==20936:
            iii=9
        for i in [0,1]: 
            # self.lengthGroups.remove(links[i])
            self.remove(links[i])
            if self.unchangedLinks.hasValue(links[i]):
                self.unchangedLinks.removeValue(links[i])
            if links[i] not in self.changedLinkInds:
                self.changedLinkInds.add(links[i])
        
        # x=self.getLinksWithWrongSize()
        # if self.shakeLinksWithWrongLength<x:
        #     iii=9
        #     raise Exception("wrong link count works incorrectly", score)
        
        point1 = tuple(sorted([score[2][0], score[2][1]])) #(A,B)
        point2 = tuple(sorted([score[2][2], score[2][3]])) #(U,V)
        self.existingLinks.remove(point1)
        self.existingLinks.remove(point2)
        point3 = tuple(sorted([score[2][1], score[2][2]])) #(B,U)
        point4 = tuple(sorted([score[2][0], score[2][3]])) #(A,V)
        self.existingLinks.add(point3)
        self.existingLinks.add(point4)
        self.links[links[0]][0], self.links[links[0]][1] =  point3[0], point3[1]
        self.links[links[1]][0], self.links[links[1]][1] =  point4[0], point4[1]

        for i in [0,1]: 
            # self.lengthGroups.add(links[i])
            self.add(links[i])

       
    def saveFile(self, fn):
        U.saveData(fn)
        print(fn, "saved")


path = "sampleData/"
fn = "data-pvalue-5-fin-min" #Original graph file to randomize !!!
print(fn, "<<")

U = UniversalDS(path + fn + ".json") #create UniversalDS instance  
print("created U")
NS = NewSmartRandomizer(U, randomizablePart=0.95, resultFn=path + fn + '-randomized' + ".json", assertions=False, chrs = U.chrs)