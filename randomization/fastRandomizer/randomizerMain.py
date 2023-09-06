from Utils import LengthGroups, DinamicList
import Utils
import random, math, copy, csv, time
import json, os
import cooler
import pandas as pd
import matplotlib.pyplot as plt

class RandLauncher:
    """
    This class reads template and launches the randomization
    """
    def __init__(self, templateName="template.json"):
        self.templateName=templateName
        self.readTemplate() #sets self.version - either "cool" or "links"
        self.randomize()
        

    def readTemplate(self):
        #reads the template file with instructions for the randomization
        templateName= self.templateName
        with open(templateName, 'r') as f:
            self.template = json.load(f)
        print("Template is read")
        if self.template["useLinkList"]:
            #link list is given
            self.version="links"
        else:
            #cool file is given
            self.version="cool"
        if not os.path.isfile(self.template["inputFn"]):
            raise Exception(f"File from template[inputFn] - {self.template['inputFn']} - not found")
        

    def randomize(self):
        #prepares all data for randomization
        if self.version=="links":
            self.randomizeLinks()
        else:
            self.randomizeCool()

    def randomizeLinks(self):
        #template asks to randomize a list of links (not .cool file)
        #first, read the links from file
        print("Reading link list..")
        with open(self.template["inputFn"], newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
            L = [[int(row[0]),int(row[1])] for row in spamreader]
        print("Link list read")
        if self.template["runNaiveRandomization"]:
            mode="randomizeLinksIgnoreLengths"
        else:
            mode=None #this means normal randomization will be launched. Used most of the time
        B = BigRandomizer(L, self.template["q"], mode=mode)
            #L - list of links. Each link is a list with 2 elements - start and end
            # q - fraction of links to swap
            # mode - usually None. If runNaiveRandomization then link lenghths will not be preserved
        
        B.saveLinks(self.template["outputFn"])
    
    def randomizeCool(self):
        
        pv = self.template["coolRelated"]["countTreshold"] #to randomize only significant links
        chr = self.template["coolRelated"]["chr"]
        inputPath = self.template["inputFn"]

        print("Reading cool file")
        cool = cooler.Cooler(inputPath)
        print("Cool file read")
        bins = cool.bins()[:].to_dict('records')
        rows = []
        pixels = cool.pixels()[:].to_dict('records')
        for px in pixels:
            if px["count"] > pv:
                a = bins[px["bin1_id"]]
                b = bins[px["bin2_id"]]
                if str(a["chrom"]) == chr and str(b["chrom"]) == chr:
                    link = ((a["start"] + a["end"]) // 2,
                            (b["start"] + b["end"]) // 2)
                    rows.append(link)
        rows = [[link[0], link[1]] for link in rows]

        if self.template["runNaiveRandomization"]:
            mode="randomizeLinksIgnoreLengths"
        else:
            mode=None #this means normal randomization will be launched. Used most of the timeprintMode=None
        B = BigRandomizer(rows, self.template["q"], mode=mode)
        B.saveLinks("./tmpRandomizedLinkList.csv") #Will be used to get cool file back

        inputPathList = "./tmpRandomizedLinkList.csv"
        inputPathCool = inputPath
        outputPath = self.template["outputFn"]

        with open(inputPathList, newline='') as f:
            reader = csv.reader(f)
            randomized = list(reader)
        cool = cooler.Cooler(inputPathCool)
        bins = cool.bins()[:].to_dict('records') #segmentu saraksts
        middle_to_id = {}
        for id, bin in enumerate(bins):
            if bin["chrom"] not in middle_to_id:
                middle_to_id[bin["chrom"]] = {}
            middle_to_id[bin["chrom"]][(bin["start"]+bin["end"])//2] = id #id in bins
        rows = []
        pixels = cool.pixels()[:].to_dict('records')
        i = 0
        existingPixels = set()
        for px in pixels:
            if px["count"] > pv:
                a = bins[px["bin1_id"]]
                b = bins[px["bin2_id"]]
                if str(a["chrom"]) == chr and str(b["chrom"]) == chr:
                    # px["bin1_id"] = middle_to_id[chr][(a["start"] + a["end"]) // 2]
                    # px["bin2_id"] = middle_to_id[chr][(b["start"] + b["end"]) // 2]
                    px["bin1_id"] = middle_to_id[chr][int(randomized[i][0])]
                    px["bin2_id"] = middle_to_id[chr][int(randomized[i][1])]
                    existingPixels.add((px["bin1_id"], px["bin2_id"], chr))
                    i += 1
        #Remove duplicates that exist due to significance - randomization might have produced links that were not significant before but were present in the cool file
        pixels = [px for px in pixels if not( (px["bin1_id"], px["bin2_id"], chr) in existingPixels and px["count"]<=pv )]

        cooler.create_cooler(outputPath, pd.DataFrame(bins), pd.DataFrame(pixels))
        print(outputPath, " saved!")


class BigRandomizer:
    def __init__(self, linkList, q, seed=None, mode=None, printMode=None):
        # linkList == [[10000000, 15400000], [], [], ...]
        # q=[0..1] percent of links to randomize 
        # seed - for random.seed
        # mode - either anything at all or "randomizeLinksIgnoreLengths" which runs a naive version
        # print mode - sets amount of information printed during randomization - either None for default or a dict with values, see self.initialization
        # self.links will contain the randomized links
        self.links = linkList
        if q==0: return
        self.linkCount = len(self.links)
        self.q = q
        self.printMode=printMode
        if seed is None:
            random.seed()
        else:
            random.seed(seed)
        print("BigRandomizer init")
        self.initialization()
        if mode is None:
            self.startSwapping()
            self.doMillenium()
        elif mode=="randomizeLinksIgnoreLengths":
            self.randomizeLinksIgnoreLengths()

    def linkLength(self,  linkInd):
        #Takes a link by index, returns link length
        #paņem linku pēc indeksa, atgriež linka garumu
        return self.links[linkInd][1] - self.links[linkInd][0]
    def normLinkLength(self,  linkInd):
        #Takes link index, returns normalized link length
        #paņem linku pēc indeksa, atgriež linka garumu
        return self.normLength(self.linkLength(linkInd))
    def normLength(self, length):
        #Takes link length, returns normalized link length
        #paņem linka garumu, atgriež linka garumu
        return length >> self.modLen
    

    def initLengthMod(self):
        #Calculates modLen - based on bin sizes. Used to approximate link lengths and deal with smaller numbers
        lenList = sorted([self.linkLength(i) for i in range(len(self.links))])
        self.modLen = Utils.bitLength(lenList[0])
        lenMap = {}
        for i in range(len(self.links)):
            v = self.normLinkLength(i)
            if v not in lenMap: 
                lenMap[v] = 0
            lenMap[v] += 1
        self.maxModInd = ((self.segments[-1] - self.segments[0]) >> self.modLen) + 1 

    def initialization(self):
        self.params = {}
        if self.printMode is None: #controls amount of information printed
            self.params["modes"] = {"show output": True,
                                    "show minimum output": True,
                                    "show graphs": True}
        else:
            self.params["modes"] = self.printMode

        self.initStats()
        self.case=0 # 0-swapping new link; 1-improving lengths

        self.originalLinks = set([(el[0],el[1]) for el in self.links]) #set of (A,B)
        self.existingLinks = set([(el[0],el[1]) for el in self.links]) #set of (A,B)

        self.segments = sorted(list(set([link[0] for link in self.links]+[link[1] for link in self.links])))
        self.segmentIndex = {self.segments[i]: i for i in range(len(self.segments))}
        minSegmentCountParam = round(math.log2(len(self.segments))) #+-32

        self.initLengthMod()
        self.lengthGroups = LengthGroups(self, minElementCount=32, minLengthCount=2)

        self.params["longestToleratedLength"] = self.lengthGroups.maxExistingLength*1.25 #TODO
        
        self.changedLinkInds = set() #set of link indeces
        unchangedLinkInds = [i for i in range(self.linkCount)]
        random.shuffle(unchangedLinkInds)
        self.unchangedLinks = DinamicList(unchangedLinkInds, 'unchangedLinks')

        self.originalLinkLengthsDraw = [link[1]-link[0] for link in self.links]

        linkLengths = [link[1]-link[0] for link in self.links]
        possibleLinkLengths = set(linkLengths)
        segments = set([link[0] for link in self.links]+[link[1] for link in self.links])
        
        #Calculating the max theoretically possible link count that can be swapped
        possibleLinksForLength = {} #garumam piekārto max iespējamo linku skaitu; 
        for linkLength in possibleLinkLengths:
            possibleLinksForLength[linkLength]=0
            for segment in segments:
                if (segment+linkLength) in segments:
                    possibleLinksForLength[linkLength]+=1
        existingLinksForLength = {}
        for le in linkLengths:
            if le not in existingLinksForLength:
                existingLinksForLength[le]=0
            existingLinksForLength[le]+=1
        blackLine = []
        for le in possibleLinkLengths:
            count = possibleLinksForLength[le] - existingLinksForLength[le]
            for i in range(count):
                blackLine.append(le)
        self.blackLine = blackLine

        linksCanRanodmize=0
        for le in possibleLinkLengths:
            linksCanRanodmize+=min((possibleLinksForLength[le] - existingLinksForLength[le]), (existingLinksForLength[le]))
        self.linksCanBeRandomized = linksCanRanodmize
        self.linksToRandomize = int(self.linksCanBeRandomized*self.q)
        print(f"Out of {self.linkCount} links, only {self.linksCanBeRandomized} can theoretically be randomized.\n\
              User asks to randomize {self.q} of links, therefore {self.linksToRandomize} links will be randomized")


    def startSwapping(self):
        #main method that starts swapping the links
        self.counter=0
        
        while len(self.changedLinkInds) < self.linksToRandomize:
            if self.lengthGroups.badLinkLengthCount<0:
                print("ER4")
            self.counter+=1
            if ((self.counter%1000) == 0):
                self.doMillenium()
                #Sometimes we choose to give up - if no link could be swapped in 5000 tries, we finish
                #If k=5 x 1000 in a row no link could be swapped, program gives up and finishes the process
                if not (hasattr(self, "swapLimitation")):
                    k=5 
                    self.swapLimitation = [0 for _ in range(k)] + [0,k]
                self.swapLimitation[self.swapLimitation[-2]]=len(self.changedLinkInds) 
                self.swapLimitation[-2]=(self.swapLimitation[-2]+1)%(self.swapLimitation[-1])
                if self.counter>1000 and len(set(self.swapLimitation[:-2]))==1: break #give up
                
            #Choice - whether improve lengths or swap new links; pasaka, vai uzlabosim, vai mainīsim
            if self.throwCoin():
                #improve lengths
                self.case=1  #0-swapping new link; 1-improving lengths
                self.swapAndImproveLinkLength()
            else:
                #swap new links
                self.case=0 #0-swapping new link; 1-improving lengths
                self.swapNewLink()
                
            self.stats["cases"][self.case]+=1 #0-swapping new link; 1-improving lengths                

    def swapAndImproveLinkLength(self):
        uncheckedLinkInd = self.lengthGroups.positiveCandidate()
        #Searching a candidate somewhere on the chromosome. meklējam kandidātu sarakstu (nekāda apkārtne)
        score = self.getBestCandidate(uncheckedLinkInd) 
        self.doTheSwap(score)

    def swapNewLink(self):
        #we need to swap a link that was never swapped
        if random.random()<0.25:
            #worsens badLinkCount
            #25% chance to choose just any unchanged link, respecting weights - biežāk ņemot linkus no grupām, kuros ir relatīvi maz samainīto linku
            uncheckedLinkInd = self.lengthGroups.getRandomElemFromRandomWeightedGroup()
        else:
            #other 50% chance to take a positive length link from a positive group
            uncheckedLinkInd = self.lengthGroups.getPositiveElemFromRandomPositiveGroup()
        score = self.getBestCandidate(uncheckedLinkInd)
        self.doTheSwap(score)
    
    def getBestCandidate(self, chosenLinkInd):
        #for a given linkInd, find a good candidate
        #by first finding several candidates and then choosing one of them
        #returns a score - list in form [A,B,U,V,chosenLinkInd,candidateLinkInd, benefitScore]
        
        allScores = [] #List of all schecked scores. We will choose one of the best afterwards. 
                        #saraksts ar visiem pārbaudītajiem scores. Pēc tam jāpaņem labākais

        if self.lengthGroups.positiveGroups.length()==0: 
            #Take 1000 ranodm links
            for i in range(1000):
                candidateLinkInd = self.lengthGroups.getRandomElemFromRandomAnyGroup()
                score = self.getScore(chosenLinkInd, candidateLinkInd)
                if not self.isUndefScore(score):
                    allScores.append(score)
        else:
            #Choose 1000 positive links
            allLinksFromPositiveGroups = self.lengthGroups.getAllLinksFromPositiveGroups()
            if len(allLinksFromPositiveGroups)>1000:
                random.shuffle(allLinksFromPositiveGroups)
                allLinksFromPositiveGroups = allLinksFromPositiveGroups[:1000]
            for candidateLinkInd in allLinksFromPositiveGroups:
                score = self.getScore(chosenLinkInd, candidateLinkInd)
                self.stats["cases pool"][self.case][score[-1]//2 + 3]+=1
                if not self.isUndefScore(score):
                    allScores.append(score)
        #At this point we got the allScores list
        #it can be empty
        if len(allScores)==0:
            return self.undefScore()
        
        allScores = sorted(allScores, key = lambda x: (x[-1]+2*x[-2]), reverse=True) #Labākais score ir priekšā.  Labākais ir tādsa, kuram summa ir labākā
        #Sort so that best scores are in the front

        #take a random; with high probability from the beginning; with lower probability with lower score
        firstHighest=len(allScores)

        index = math.floor(((random.random())**4)*(firstHighest)) #random.random can not generate 1.0; it generates [0.0...1.0)], hence no -1
        #index - chosen in the sorted allScores
        bestScore = allScores[index]

        self.stats["case scores"][self.case][bestScore[-1]//2 + 3]+=1
        return bestScore

    def getScore(self, chosenLinkInd, candidateLinkInd):
        #chosenLinkInd, candidateLinkInd are link indeces. We need to calculate a swapping score for this swap
        
        a = self.links[chosenLinkInd]
        b = self.links[candidateLinkInd]
        #First check if they have a common endpoint
        if a[0]==b[0] or a[0]==b[1] or a[1]==b[0] or a[1]==b[1]:
            #bad case
            return self.undefScore()
        
        #a[0]-b[1] un a[1]-b[0] ir jaunie linki - new links
        score1 = self.getScoreForSwap([a[0], a[1], b[0], b[1], chosenLinkInd, candidateLinkInd])
        
        #a[0]-b[0] un a[1]-b[1] ir jaunie linki - new links
        score2 = self.getScoreForSwap([a[0], a[1], b[1], b[0], chosenLinkInd, candidateLinkInd])

        #return a better score
        return score1 if score1[-1]>score2[-1] else score2
    
    def getScoreForSwap(self, p):
        # p == [A,B,U,V,chosenLinkInd,candidateLinkInd];  which creates links AV and BU
        #returns [A,B,U,V,chosenLinkInd,candidateLinkInd, WEIGHT, benefitScore] where positive benefitScore is better  

        newPoints = [tuple(sorted([p[0], p[3]])), tuple(sorted([p[1], p[2]]))] #[(A,V), (B,U)]
        if newPoints[0] in self.existingLinks or newPoints[1] in self.existingLinks:
            #creating a link that already exists
            return self.undefScore()
        if newPoints[0] in self.originalLinks or newPoints[1] in self.originalLinks:
            #creating a link that was in the original graph
            return self.undefScore()
        if self.normLength(newPoints[0][1]-newPoints[0][0])>self.params["longestToleratedLength"] or\
            self.normLength(newPoints[1][1]-newPoints[1][0])>self.params["longestToleratedLength"]:
            #link is too long
            return self.undefScore()
        
        oldLenGroup1 = self.lengthGroups.getGroupForLinkInd(p[4])
        oldLenGroup2 = self.lengthGroups.getGroupForLinkInd(p[5])

        newLength1 = self.normLength(newPoints[0][1]-newPoints[0][0])
        newLength2 = self.normLength(newPoints[1][1]-newPoints[1][0])
        newLenGroup1 = self.lengthGroups.getGroup(newLength1)
        newLenGroup2 = self.lengthGroups.getGroup(newLength2)

        benefitScore = 0

        groupDif = {}
        for groupInd in [oldLenGroup1, oldLenGroup2, newLenGroup1, newLenGroup2]:
            if groupInd not in groupDif:
                groupDif[groupInd]=self.lengthGroups.getGroupDifference(groupInd)
        
        #calculate benefit score. If it is positive, we improve 
        osign = lambda x: 1 if x > 0 else -1 #old sign. 
        for groupInd in [oldLenGroup1, oldLenGroup2]:
            benefitScore += osign(groupDif[groupInd])
            groupDif[groupInd]-=1
        
        nsign = lambda x: 1 if x < 0 else -1 #signfunction for new lengths
        for groupInd in [newLenGroup1, newLenGroup2]:
            benefitScore += nsign(groupDif[groupInd])
            groupDif[groupInd]+=1
        
        #At this point we have calculated benefitScore
        #Now we calculate bonusScore, a coeficient that promotes swapping of short links (links from groups with many links)
        bonusForCandidateLink = self.getWeight(oldLenGroup2)
        bonusForNewLink1 = self.getWeight(newLenGroup1)
        bonusForNewLink2 = self.getWeight(newLenGroup2)
        totalBonus = bonusForCandidateLink + bonusForNewLink1 + bonusForNewLink2 #in range 0..3. we want 3,2,1,0
        
        p.append(totalBonus)
        p.append(benefitScore)
        return p
    
    def getWeight(self, groupInd):
        #weight that promotes choosing shorter links for a swap
        #oldLenGroup2 = self.lengthGroups.getGroupForLinkInd(linkInd)
        # #Now we calculate bonusScore, a coeficient that promotes swapping of short links (links from groups with many links)
        # #density = (linku skaits ar doto garumu)                    / (linka garums (normētais))
        # density = (self.lengthGroups.groupInitLength[oldLenGroup2])/(self.normLinkLength(linkInd))
        # #ratio = (cik neesam samainījuši no linkiem grupā)        /(linki grupā)
        # ratio = (self.lengthGroups.firstSwappedLink[oldLenGroup2])/(self.lengthGroups.groupInitLength[oldLenGroup2])
        # weight = density*ratio

        #weight = 1 - [stabiņa daļa zem zaļās līknes]
        #stabiņa daļa zem zaļās līknes = [samainīto linku skaits grupā] / [linku skaits grupā]
        if self.lengthGroups.groupInitLength[groupInd]==0: return 0
        weight = self.lengthGroups.firstSwappedLink[groupInd] / self.lengthGroups.groupInitLength[groupInd] #kāda daļa ir nemainītie - fraction of links unchanged
        return weight*weight
        
    ############################
    def startImproving(self):
        #after enough links were swapped, we start to improve link length distribution
        print("\nSTART IMPROVING LINK LENGTHS")
        minImprovementCount = 5 #out of 1000
        self.counter=1
        currentBadLinks = self.lengthGroups.badLinkLengthCount
        #the (self.lengthGroups.badLinkLengthCount - currentBadLinks) > minImprovementCount is checked on each 1000th iteration
        while ((self.counter%1000)!=0 or (currentBadLinks-self.lengthGroups.badLinkLengthCount) > minImprovementCount):
            self.swapAndImproveLinkLength()
            self.case=1
            self.stats["cases"][self.case]+=1 #0-swapping new link; 1-improving lengths
            if (self.counter%1000)==0:
                currentBadLinks = self.lengthGroups.badLinkLengthCount
                self.doMillenium()
            self.counter+=1
        
        #get rid of the longest links
        positiveLinkInd = self.lengthGroups.positiveCandidate()

    def throwCoin(self):
        #0 or false - swap new link
        # We want this to happen at least 75% of time but no more than 90%        
        p = 1 - (len(self.changedLinkInds)/(self.linksToRandomize)) #0..1 (actually -inf to 1)
        p = 0.75 + 0.15*p
        coin = random.random()
        if coin>p:
            return True
            #return "improvment of lengths" #happens 10-25% of the time
        else:
            return False
            #return "swapping new links" #will happen often at first, and less frequently later

    def isUndefScore(self, score):
        return score[-1] == -6

    def undefScore(self):
        #[A,B,U,V,ind1,ind2,weight,benefitScore]
        return [-1, -1, -1, -1, -1, -1, -1, -6]

    def doTheSwap(self, score):
        #Performs the swap of two links and maintains all other lists
        # score is [A,B,U,V,chosenLinkInd,candidateLinkInd, benefitScore] for swap AV and BU
        if self.isUndefScore(score): return 0        
        linkInds = score[4:6]  #link indeces
        for i in [0,1]: 
            linkInd = linkInds[i]
            # self.remove(linkInd)
            self.lengthGroups.remove(linkInd)
            if self.unchangedLinks.hasValue(linkInd):
                self.unchangedLinks.removeValue(linkInd)
                self.changedLinkInds.add(linkInd)

        point1 = tuple(sorted([score[0], score[1]])) #(A,B)
        point2 = tuple(sorted([score[2], score[3]])) #(U,V)
        self.existingLinks.remove(point1)
        self.existingLinks.remove(point2)
        point3 = tuple(sorted([score[1], score[2]])) #(B,U)
        point4 = tuple(sorted([score[0], score[3]])) #(A,V)
        self.existingLinks.add(point3)
        self.existingLinks.add(point4)
        self.links[linkInds[0]][0], self.links[linkInds[0]][1] =  point3[0], point3[1]
        self.links[linkInds[1]][0], self.links[linkInds[1]][1] =  point4[0], point4[1]
        for i in [0,1]: 
            self.lengthGroups.add(linkInds[i])

    def testSwapLinks(self):
        #Data integrity check, used in development
        self.allNeverSwappedLinkInds = set([i for i in range(len(self.links))])
        def dataIntegrityCheck():
            for link in self.links:
                if not (link[0]<link[1]): 
                    print("ER1")
                #linkGroup = vai elInd no linka ir mazāks vai vienāds par group elem garumu. Varbūt ielikts nepareizā grupā
                #Vai group elem no group no elemInd no link == linkInd
                linkInd = self.links.index(link)
                
                trueGroup = self.lengthGroups.getGroupForLinkInd(linkInd)
                groupBy = self.lengthGroups.elemGroup[linkInd]
                if not (trueGroup==groupBy): 
                    print("ER2")
            for groupID in self.lengthGroups.groupElem.keys():
                group = self.lengthGroups.groupElem[groupID] #a list
                firstSwappedLink = self.lengthGroups.firstSwappedLink[groupID]
                onlyHadUnswapped=True
                for i in range(len(group)):
                    if (i<firstSwappedLink) == (group[i] in self.allNeverSwappedLinkInds):
                        pass
                    else:
                        print("ER3")

        dataIntegrityCheck() 

        for i in range(1000):
            linkInd1, linkInd2 = random.randint(0, len(self.links)-1), random.randint(0, len(self.links)-1)
            self.lengthGroups.remove(linkInd1)
            self.lengthGroups.remove(linkInd2)
            self.links[linkInd1][1], self.links[linkInd2][1] = self.links[linkInd2][1], self.links[linkInd1][1]
            for linkInd in [linkInd1, linkInd2]:
                if self.links[linkInd][1]<self.links[linkInd][0]: self.links[linkInd] = [self.links[linkInd][1],self.links[linkInd][0]]
            self.allNeverSwappedLinkInds.discard(linkInd1)
            self.allNeverSwappedLinkInds.discard(linkInd2)
            self.lengthGroups.add(linkInd1)
            self.lengthGroups.add(linkInd2)

        dataIntegrityCheck() 

    def plotRndProgress(self, originalLinkLengthsDraw, currentLinkLengthsDraw, swappedLinkLengthsDraw, n):
        def plotLengths(ax, originalLinkLengthsDraw, currentLinkLengthsDraw, swappedLinkLengthsDraw):
            #MAX_DRAWABLE_LINK=0.5*1e8
            MAX_DRAWABLE_LINK = max(currentLinkLengthsDraw)
            ax.hist([l for l in originalLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=int(MAX_DRAWABLE_LINK//10e5), alpha=0.9, label="Original", histtype='step')
            ax.hist([l for l in currentLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=int(MAX_DRAWABLE_LINK//10e5), alpha=0.9, label="Randomized", histtype='step')
            ax.hist([l for l in swappedLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=int(MAX_DRAWABLE_LINK//10e5), alpha=0.9, label="Swapped", histtype='step')
            # ax.hist([l for l in originalLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=list(sorted(list(set(originalLinkLengthsDraw)))), alpha=0.9, label="Original", histtype='step')
            # ax.hist([l for l in currentLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=list(sorted(list(set(currentLinkLengthsDraw)))), alpha=0.9, label="Randomized", histtype='step')
            # ax.hist([l for l in swappedLinkLengthsDraw if l<=MAX_DRAWABLE_LINK], bins=list(sorted(list(set(swappedLinkLengthsDraw)))), alpha=0.9, label="Swapped", histtype='step')
            # ax.hist([l for l in self.blackLine if l<=MAX_DRAWABLE_LINK], bins=list(sorted(list(set(self.blackLine)))), alpha=0.9, label="extra available", histtype='step', color="black")

            #add black line that show how many links of a given length can be added to a graph till it is full


        def plotProgress(ax, totalLinkCount, randomizableLinkCount, linksRandomized, badLinkCount):
            names = ["Links total", "Links to rnd", "Links rnd-ed", "Bad links"]
            counts = [totalLinkCount, randomizableLinkCount, linksRandomized, badLinkCount]
            ax.bar(names, counts, label=names, color=["forestgreen","limegreen", "turquoise", "teal" ])
            y_pos = range(len(names))
            ax.set_xticks(y_pos, names, rotation=90)
            #ax.margins(y=0.2)
        
        # originalLinkLengthsDraw = lengths
        # currentLinkLengthsDraw = [l*2 for l in lengths]
        fig, ax = plt.subplot_mosaic("AAAB;AAAB;AAAB")
        #plt.tight_layout(pad=4.0)
        ax["A"].set_xlabel("Link length")
        ax["A"].set_ylabel("Count of links")
        plotLengths(ax["A"], originalLinkLengthsDraw, currentLinkLengthsDraw, swappedLinkLengthsDraw)
        ax["A"].legend()
        plotProgress(ax=ax["B"], totalLinkCount=len(originalLinkLengthsDraw),\
                    randomizableLinkCount=self.linksToRandomize,\
                    linksRandomized=len(self.changedLinkInds),\
                    badLinkCount=self.lengthGroups.badLinkLengthCount)
        fig.suptitle(f"After {n} iteration")
        plt.tight_layout(pad=1.5)
        path = "rez"
        isExist = os.path.exists(path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path)
        plt.savefig(f"rez/rnd_progress_after_{n}.png", dpi=800)
        plt.clf()
        plt.close()
    def initStats(self):
        self.stats = {
            "iterations": 0,
            "time": 0,
            "cases": [0,0],
            "positive groups": 0,
            "positive links": [0,0], #["never swapped links in positive groups; swapped links in positive groups"]
            "bad links": 0,
            "swapped links": 0,
            "case scores": [[0,0,0,0,0,0],[0,0,0,0,0,0]],
            "cases pool": [[0,0,0,0,0,0],[0,0,0,0,0,0]],
        }
        self.lastCheckpointTime = time.time()
    def updateStats(self):
        self.stats["iterations"] = self.counter 
        curTime = time.time() 
        self.stats["time"] = round(curTime - self.lastCheckpointTime, 2)
        self.lastCheckpointTime = curTime    
        self.stats["positive groups"]= len(self.lengthGroups.positiveGroups.data)
        self.stats["bad links"] = self.lengthGroups.badLinkLengthCount
        self.stats["swapped links"] = len(self.changedLinkInds)     
        
        for positiveGroupInd in self.lengthGroups.positiveGroups.data: #["never swapped links in positive groups; swapped links in positive groups"]
            group = self.lengthGroups.groupElem[positiveGroupInd]
            self.stats["positive links"][0]+=self.lengthGroups.firstSwappedLink[positiveGroupInd]
            self.stats["positive links"][1]+=len(group)-self.lengthGroups.firstSwappedLink[positiveGroupInd]

    def doMillenium(self):
        #after 1000 iters, print everything
        if self.params["modes"]["show output"]:
            self.updateStats()
            if self.params["modes"]["show minimum output"]:
                print(f'Iterations: {self.stats["iterations"]};  time:{self.stats["time"]}; already swapped: {self.stats["swapped links"]} / {self.linksToRandomize};')
            else:
                print(self.stats)
            self.dumpStats()
            if self.counter%20000 ==0:
                self.saveStatsToFile()
            self.initStats()
        
        if self.params["modes"]["show graphs"]:
            if ( self.counter>0 and self.counter%10000)==0:
                originalLinkLengthsDraw = self.originalLinkLengthsDraw
                currentLinkLengthsDraw = [link[1]-link[0] for link in self.links]
                swappedLinkLengthsDraw = [self.links[i][1] - self.links[i][0] for i in range(self.linkCount) if i in self.changedLinkInds]
                n = self.counter
                self.plotRndProgress(originalLinkLengthsDraw, currentLinkLengthsDraw, swappedLinkLengthsDraw, n)

    def dumpStats(self):
        #dumps one row to stats file
        header = []
        row = []
        for key in self.stats.keys():
            if type(self.stats[key])==type(100) or type(self.stats[key])==type(100.1) or type(self.stats[key])==type("st"):
                header.append(key)
                row.append(self.stats[key])
            elif type(self.stats[key])==type({"k":12}):
                for k in self.stats[key].keys():
                    header.append(k)
                    row.append(self.stats[key][k])
            elif type(self.stats[key])==type([]): #list or list of lists
                i=-1
                for el in self.stats[key]:
                    i+=1
                    if type(el)==type([]): #list of lists
                        j=-1
                        for e in el:
                            j+=1
                            header.append(f"{key}-{i}{j}")
                            row.append(self.stats[key][i][j])
                    else:
                        #list of int
                        header.append(f"{key}-{i}")
                        row.append(self.stats[key][i])
        
        if not (hasattr(self, "statsDump")):
            self.statsDump = []
            self.statsDump.append(header)
        self.statsDump.append(row)
    
    def saveStatsToFile(self, fn="./dumpedStats.csv"):
        return
        with open(fn, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            spamwriter.writerows(self.statsDump)
        print (fn, "saved")
    
    def saveLinks(self, fn="./dumpedLinks.csv"):
        with open(fn, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            spamwriter.writerows(self.links)
        print (fn, "saved")

    def randomizeLinksIgnoreLengths(self):
        #For result validation, the naive version of the randomization
        self.existingLinks = set([tuple(sorted(link)) for link in self.links])
        while (len(self.existingLinks) - self.linkCount) < self.linksToRandomize:
            linkInds = [random.randint(0, self.linkCount-1), random.randint(0, self.linkCount-1)]
            if len(set([self.links[linkInds[0]][0], self.links[linkInds[1]][0], self.links[linkInds[0]][1], self.links[linkInds[1]][1]])) == 4:
                index = [random.randint(0, 1), random.randint(0, 1)]
                link = [copy.copy(self.links[linkInds[0]]), copy.copy(self.links[linkInds[1]])]
                link[0][index[0]], link[1][index[1]] = link[1][index[1]], link[0][index[0]]
                link[0] = sorted(link[0])
                link[1] = sorted(link[1])
                if tuple(link[0]) not in self.existingLinks and tuple(link[1]) not in self.existingLinks:
                    self.links[linkInds[0]] = link[0]
                    self.links[linkInds[1]] = link[1]
                    self.existingLinks.add(tuple(link[0]))
                    self.existingLinks.add(tuple(link[1]))
            if len(set([tuple(link) for link in self.links]))!=len(self.links):
                iii=9
        iii=9