# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 07:40:10 2023

@author: asizo
"""
from universal import UniversalDS, ChrData
#import universal 
import csv
import json

 
class Cliques:
    def __init__(self, owner,  minC3TissueCount=1, tissueMask=0):
        #calculates cliques - triangles of length 3 and applies filters to triangles
        self.owner = owner #ChrData instance
        self.DS = owner.DS
        self.minC3TissueCount = minC3TissueCount
        self.tissueMask = tissueMask
        self.ch = self.owner.ch
        if type(tissueMask)==str: 
            try:
                self.tissueMask = self.owner.owner.tissueBits[tissueMask]
                #maybe tissue name was  instead of bits
            except KeyError:
                self.owner.error("tissueMask is invalid", tissueMask)
        elif type(tissueMask)==list:
            try:
                m=0
                for tisName in tissueMask:
                    m = (m | (self.owner.owner.tissueBits[tisName]))
                self.tissueMask = m
                #maybe list of tissue names was given
            except KeyError:
                self.owner.error("tissueMask is invalid", tissueMask)
        self.getTriangles(self.minC3TissueCount, self.tissueMask)
    
    def getTriangles(self, minC3TissueCount=1, tissueMask=0, owner=None):
        #sets and gets list of all triangles that pass the filters
        def getBitCount(bits): #returns number of 1s in binary form
            count = 0
            while (bits):
                bits &= (bits-1)
                count+= 1
            return count

        if owner is not None: self.owner = owner #owner changes
        allCliques = self.owner.getC3() #list of tuples (A, B, C, bit)
        
        allCliques = [cl for cl in allCliques if getBitCount(cl[-1])>=self.minC3TissueCount and ((cl[-1]&self.tissueMask)==self.tissueMask)]
        self.allCliques = allCliques
        
        return self.allCliques
    
    def getAllSegments(self, cliques=None):
        #from list of cliques get list of segments, by indeces
        if cliques is None: cliques = self.allCliques
        
        s = set()
        for cl in cliques:
            for i in [0,1,2]:
                s.add(cl[i])
        return sorted(list(s))
    
    def getAllSegmentsLoci(self):
        #from instance of cliques get list of segments, by loci
        indeces = self.getAllSegments()
        loci = self.owner.listOfIndsToListOfLoci(indeces)
        return loci
    
    def getListOfSegmentLociSegments(self):
        #returns list of nodes, each node in form [A, B] where A and B are loci of the segment
        segInds = self.getAllSegments()
        return [self.owner.segments[i] for i in segInds]
    
    def getLinks(self):
        #returns list of (A,B) with A<B and A-B is a link in at least 1 calculated CL3
        #uses self.allCliques
        links = set()
        for cl3 in self.allCliques:
            #cl3 is in form (A,B,C,bit)
            links.add((cl3[0], cl3[1]))
            links.add((cl3[0], cl3[2]))
            links.add((cl3[1], cl3[2]))
        return sorted(list(links), key=lambda x: (int(x[0]), int(x[1])) )
    
    def getLinksList(self):
        #returns list of (A,B) with A<B and A-B is a link in at least 1 calculated CL3. With repetitions
        #uses self.allCliques
        links = []
        for cl3 in self.allCliques:
            #cl3 is in form (A,B,C,bit)
            links.append((cl3[0], cl3[1]))
            links.append((cl3[0], cl3[2]))
            links.append((cl3[1], cl3[2]))
        return sorted(links, key=lambda x: (int(x[0]), int(x[1])) )
    
    def getC3Counts(self):
        #returns {"tis1": 1234, "tis2": 2345, "tis3":C3CountInTissue..} 
        #self.allCliques is a list of (A,B,C,bit)
        U = self.owner.owner
        rez = {tis:0 for tis in U.tissues}
        for cl3 in self.allCliques:
            tisList = U.getTissueNameList(cl3[3])
            for tis in tisList:
                rez[tis]+=1
        return rez
        
        
class Bases:
    def __init__(self, owner, withOr=True):
        #oneChrLinks == list of [Aind, Bind, bits]
        self.withOr = withOr #true if making "OR-bases", false if making "AND-bases". It changes the way of counting base degree slightly
        self.owner = owner #ChrData object/ instance
        self.DS = owner.DS
        self.N=2
        self.links = self.owner.links
        self.ch = self.owner.ch
        self.tissueBits = self.owner.owner.tissueBits #{"aCD4": 1, "EP": 2, "Ery": 4, "FoeT": 8,...
        self.cliques = self.owner.getC3() # set of (A, B, C, bits)
        self.adjCandidates = self.calculateBaseAdj()
        self.allBases = self.setBases() #sets self.bases this can be commented out - each function calls this again to have most up-to-date data
        
    
    def calculateBaseAdj(self):
        # creates adj of bases. A base is a link that is a part of many triangles. Alt. many triangles "sit" on a base
        # (U, Z) : { (V, C3tisbit(V)), (), ... }
        #C3tisbit(V) is a bitmap of those tissues, that are present in all 3 triangle links
        
        #cliques are stored in self.cliques
        adjCliqueOfSets = dict() # keys are links (Aind, Bind)
        adjCliqueOfLists = dict()
        for clique in self.cliques:
            #clique==(A,B,C,bit)                
            for pair in [(0,1,2), (0,2,1), (1,2,0)]:  #(pair, third element) - base can be AB, BC or AC
                link = (clique[pair[0]], clique[pair[1]])
                if link not in adjCliqueOfSets: adjCliqueOfSets[link] = set()
                adjCliqueOfSets[link].add((clique[pair[2]], clique[-1]))
                  
                #this time with lists and not sets. Results should be the same
                if link not in adjCliqueOfLists: adjCliqueOfLists[link] = []
                adjCliqueOfLists[link].append((clique[pair[2]], clique[-1]))
                #got adj structure for bases - for each link (U,Z) got adjacent nodes (U, Z) : { (V, C3tisbit(V)), (), ... }
        
        if self.owner.owner.assertion:
            print("Asserting validity of calculateBaseAdj results")
            #check sakrit garumi
            if len(adjCliqueOfSets)!=len(adjCliqueOfLists) : 
                self.owner.owner.warning("unequal lengths in Bases.calculateBaseAdj when calculating in 2 different ways", [adjCliqueOfSets, adjCliqueOfLists])
            for key in adjCliqueOfLists.keys():
                if set(adjCliqueOfLists[key])!=adjCliqueOfSets[key]: 
                    self.owner.owner.error("unequal sets with lists Bases.calculateBaseAdj", [adjCliqueOfSets, adjCliqueOfLists])
        #baseLinks = [link for link in adjCliqueOfSets.keys() if len(adjCliqueOfSets[link]) >= 3]
        return adjCliqueOfLists #candidates
    
    def setBases(self):
        #returns dict of (A,B): (N, tissueBits), where N is the degree of base (A,B) and tissueBits is the bitmap (either OR or AND version) of 
        bases = {tis: [] for tis in self.tissueBits.keys()}
        # tissueBits - sajos audu tipos ir bazes links
        for link in self.adjCandidates.keys():
            listOfTuples = self.adjCandidates[link] # { (V, C3tisbit(V)), (), ... } for (U, Z)
            N = len(listOfTuples)
            if N<2: continue #skip bases1
            if self.owner.owner.assertion:
                if len(listOfTuples)!=len(set([v[0] for v in listOfTuples])):
                    self.owner.owner.error("One segment of Adj has different tissues in different places, should be impossible", [link, listOfTuples])
            for tisName in self.tissueBits.keys(): #"tis1", "tis2"
                tisBit = self.tissueBits[tisName]
                n = len([v[0] for v in listOfTuples if (v[1]&tisBit)>0])
                # current link is a part of n triangles in the current tissue type. There are a total of N different triangles containing current link across all tissue types
                if self.withOr:
                    if n>0:
                        bases[tisName].append([link[0], link[1], N])
                elif n>=2:
                    bases[tisName].append([link[0], link[1], n]) #and version. Note the lowercase vs uppercase N and n

        self.bases = bases
        return bases
    
    def setBasesOr(self): 
        #set and return OR-bases. Degrees are generally bigger. For any one link, base degree is the same in any tissue
        # OR-base - takes all triangles in any tissue type and count number of triangles for each link (possible base)
        # tissue bit for the base - at least one triangle is present in each of the tissues from tissue bitmap
        self.withOr=True
        self.setBases()
        return self.bases
    
    def setBasesAnd(self):
        #set and return AND-bases. Degrees are generally smaller. For any one link, base degrees are different for different tissues
        # AND-base (realisation is different, but result is the same) - take a graph with only links in one tissue type, find bases there for this tissue only
        self.withOr=False
        self.setBases()
        return self.bases
                    
    def getDegreeNBases(self, N):
        #For all tissues: using self.bases, find those bases that have degree >=N
        self.N = N
        return { tis: [triple for triple in self.bases[tis] if triple[2]>=N] for tis in self.bases.keys()  } 
    
    def getSegmentStructureInds(self, N, withOr=None, filename=None): 
        # formerly getSandraStructure
        # for every tissue, return list of segment indeces that are found to be bases of degree N+
        # if filename is set, save there
        #filename==noFileSave to not save the file
        if withOr is None: withOr = self.withOr
        if type(withOr)!=type(True): self.owner.error("withOr is not a bool", withOr)
        self.setBases() #sets self.bases
        bases = self.getDegreeNBases(N) #only those bases that are a part of at least N triangles
        rez = dict.fromkeys(bases.keys())
        for tis in bases.keys():
            segs = set()
            for triple in bases[tis]: #triple is [A, B, degree] - base (A,B) and its degree
                for i in [0,1]:
                    segs.add(triple[i])
            rez[tis] = sorted(list(segs))
        
        #saving the result and returning json
        if filename=="noFileSave": return rez
        if filename is None: filename = "baseSegmentIndsStructure.json"
        if filename.find(".json")<0: filename=filename+".json"
        out_file = open(filename, "w")
        json.dump(rez, out_file)
        out_file.close()  
        return rez
    
    def getSegmentStructureLoci(self, N, withOr=True, filename=None):
        #writes all segment loci (not indeces) to file that are part of any base with degree N+
        r = self.getSegmentStructureInds(N, withOr, "noFileSave")
        newR = {t: [] for t in r.keys()}
        for tis in r.keys():
            newR[tis] = [self.owner.segmentIndToMidpoint[el] for el in r[tis]]
        #translated segment indeces to loci. Now returning result and saving to file
        
        if filename=="noFileSave": return newR
        if filename is None: filename = "baseSegmentLociStructure.json"
        if filename.find(".json")<0: filename=filename+".json"
        out_file = open(filename, "w")
        json.dump(newR, out_file)
        out_file.close()  
        return newR
    
    def getAllSegments(self, N=None, withOr=True, tissue=None):
        #either get list of segments for bases degree N+ for all tissues if tissue is not set
        #..or for one particular tissue
        s = set()
        if N is None: N=self.N
        if N is None: N=2
        if tissue is None: #got no tissues - calculating segments for all tissues
            r = self.getSegmentStructureInds(N, withOr, "noFileSave")
            for key in r.keys():
                s.update(r[key])
        else:
            if type(tissue) == type([]): tissues = tissue #got a list of tissues
            else: tissues = [tissue] #got one tissue name
            for tis in tissues:
                if tissue in r:
                    s.update(r[tissue])
                else:
                    self.owner.error("Tissue not found in call Bases.getAllSegments", [tissue, r.keys()])
        return sorted(list(s))
    def getAllSegmentsLoci(self):
        #returns list of all segments in form of loci
        segInds = self.getAllSegments()
        return self.owner.listOfIndsToListOfLoci(segInds)
    
    def getListOfSegmentLociSegments(self):
        #returns list of nodes, each node in form [A, B] where A and B are loci of the segment
        segInds = self.getAllSegments()
        return [self.owner.segments[i] for i in segInds]
 

    
                    


class BasesOfBases:
    def __init__(self, owner, links = None):
        # owner - ChrData instance
        #links - ChrData.links by default or something else
        #only OR-bases are calculated, tissue type is ignored
        if links is None: links = owner.links
        self.links = links
        self.owner = owner
        self.DS = owner.DS
        self.chData = ChrData(owner.owner, owner.ch, self.links) #copy of a passed ChrData object instance (owner)
        self.withOr = True
        self.ch = self.owner.ch
        
    
    def reduce(self, deg):
        #calculates bases for self.links. 
        # keeps only those links that are part of N+ triangles (are bases with degree N+). N:=deg
        # modifies class instance data
        # If called twice in a row, bases of bases and basesOfBasesOfBases will be calculated
        Bor = Bases(self.chData, withOr=self.withOr)
        baseLinks = [v[:3] for v in self.links if tuple([v[0], v[1]]) in Bor.adjCandidates.keys()]
        baseDegrees = [len(Bor.adjCandidates[tuple([v[0], v[1]])]) for v in baseLinks]
        
        newLinks = []
        for i in range(len(baseLinks)):
            if baseDegrees[i]>=deg: newLinks.append(baseLinks[i])
        # newLinks now is a list of [A,B,bitmap] of only those links that are bases N
        
        self.baseNCount = len(newLinks) 
        print ("Bases {deg} of bases calculated. {n} bases found".format(deg=deg, n=self.baseNCount))
        self.links = newLinks
        self.chData.links = self.links
        
        return self.links
    
    def calc(self, deg):
        #calculates bases but does not modify instance data
        Bor = Bases(self.chData, withOr=self.withOr)
        baseLinks = [v[:3] for v in self.links if tuple([v[0], v[1]]) in Bor.adjCandidates.keys()]
        baseDegrees = [len(Bor.adjCandidates[tuple([v[0], v[1]])]) for v in baseLinks]
        
        newLinks = []
        for i in range(len(baseLinks)):
            if baseDegrees[i]>=deg: newLinks.append(baseLinks[i])
        
        count = len(newLinks)
        self.curCount = count #how many bases with degree deg are found
        print ("Bases {deg} of bases calculated without modifying the instance. {n} bases found".format(deg=deg, n=count))
        
        return newLinks
    
    def getSegments(self, links=None):
        # gets current list of segments in calculated bases
        if links is None: links = self.links
        segs = set()
        try:
            for el in links:
                segs.add(el[0])
                segs.add(el[1])
            segs = sorted(list(segs))
        except IndexError:
            self.owner.owner.error("In BasesOfBases.getSegments link is not in form [A,B,bitmap]", el)
        return segs
 
    def getAllSegmentsLoci(self):
        #returns list of all segments in form of loci
        segInds = self.getSegments()
        return self.owner.listOfIndsToListOfLoci(segInds)
    
    def getListOfSegmentLociSegments(self):
        #returns list of nodes, each node in form [A, B] where A and B are loci of the segment
        segInds = self.getSegments()
        return [self.owner.segments[i] for i in segInds]
        
    

