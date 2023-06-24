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
        self.fn = fn
        self.readData() #sets self.data, self.{chrs, tissueBits, bitToTisname}
        self.assertion = assertion #False by default. If true, validates that received data is of correct format
    
    def readData(self):
        f = open(self.fn)
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
        
        
        
        
        
    

