#!/usr/bin/env python 
# -*- coding: utf-8 -*- import sys
import math, random, json, datetime, sys, os, copy

def warningMSG(text, value):
    print("### WARNING ###", text, "###", value)
def errorMSG(text, value):
    print("\n### ERROR ###", text, "###", value)
    exit()
def existFile(PATH): return (os.path.isfile(PATH) and os.access(PATH, os.R_OK))
def existDir(PATH): return os.path.exists(PATH)
def createDir(path):
    # Check whether the specified path exists or not
    if not existDir(path):
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("### The new directory is created:", path)
def createFile(path):
    f = open(path, 'a+')
    print('### The new file is created:', path)
    f.close()
def getFileExt(path, ext):
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        errorMSG('Path do not exist', path)
    return [v for v in os.listdir(path) if v.find(ext) >= 0]
def readJson(json_file):
	json_data = open(json_file)
	values = json.load(json_data)
	json_data.close()
	print ("DONE load ", json_file)
	return values

def int2Strig(i):
    if i < 10: return "0" + str(i)
    return str(i)
# atgriez datumu no datetime
def toDMY(s):
    if isinstance(s, datetime.datetime):  return s.strftime('%d.%m.%Y.')
    return s
def date(dt): return int(dt.strftime("%d"))
# atgriez indeksu list, kur ir pirma vertiba value, -1 ja vispar nav
def indexByValue(list, value): return next((i for i in range(len(list)) if list[i] == value), -1)
def indexByTagValue(list, tag, value): return next((i for i in range(len(list)) if tag in list[i] and list[i][tag] == value), -1)
def listOfDicToTagDic(values, tag):
    dic = {}
    for v in values:
        dic[v[tag]] = v
    return dic
def listOfDicToTagPairDic(values, tags):
    dic = {}
    for v in values:
        dic[v[tags[0]] + ' ' + v[tags[1]]] = v
    return dic
def rowToDic(header, values):
    dic = {}
    for i in range(len(header)): dic[header[i]] = values[i]
    return dic
def swapList(list):
    n = len(list)
    for i in range(n):
        k = random.randint(0, n - 1)
        list[k], list[0] = list[0], list[k]
    return list
def round025(v):
    w = round(v, 2)
    if v * 4 == int(v * 4): return w
    return round(v, 1)
def index_2d(data, search):
    for i, e in enumerate(data):
        try:
            return i, e.index(search)
        except ValueError:
            pass
    return [-1, -1]
def bitLength(x):
    return 0 if x < 2 else (x + 1).bit_length()

class DinamicList:
    # Each element of data supposed to be uniq.
    def __init__(self, data, txt):
        self.data = copy.copy(data)
        self.textMSG = txt
        self.dataMap = {data[i]: i for i in range(len(data))}
        if len(self.data) != len(self.dataMap.keys()):
            warningMSG(self.textMSG + '::init:: Dinamic list data has same values', [])
            return  
    def __str__(self) -> str:
        return "data: {}".format(str(self.data))      
    def hasValue(self, value):
        return value in self.dataMap
    def randValue(self):
        if len(self.dataMap) == 0:
            warningMSG(self.textMSG + '::Dinamic list is empty', [])        
            return
        return self.data[random.randrange(0, len(self.data))]
    def addValue(self, value):
        if self.hasValue(value):
            warningMSG(self.textMSG + '::addValue:: Dinamic list try add same values', [value, self.data])        
            return
        self.dataMap[value] = len(self.data)
        self.data.append(value)
    def removeValue(self, value):
        if not self.hasValue(value):
            warningMSG(self.textMSG + '::removeValue:: Dinamic list has not value', [value, self.data])
            return
        lastValue = self.data.pop()
        if lastValue != value:
            self.data[self.dataMap[value]] = lastValue
            self.dataMap[lastValue] = self.dataMap[value]
        self.dataMap.pop(value)
    def length(self):
        return len(self.data)
    def isValid(self):
        valid = True
        if len(self.data) != len(self.dataMap.keys()):
            warningMSG(self.textMSG + '::isValid:: Dinamic list data has same values', [])
            valid = False
        n = len(self.data)
        i = 0
        while i < n and valid:
            value = self.data[i]
            if self.dataMap[value] < 0 or self.dataMap[value] >= n or self.dataMap[value] != i: 
                warningMSG(self.textMSG + '::isValid:: Dinamic list wrong dataMap', [i, value, self.dataMap[value]])
                valid = False
                break
            i += 1
        return valid

class LengthGroups:
    def __init__(self, owner, minElementCount=16, minLengthCount=4, groupWidth=200000, txt='Unknown'):
        self.owner = owner 
        self.minElementCount = minElementCount
        self.minLengthCount = minLengthCount
        self.groupWidth = self.owner.normLength(groupWidth) #If groups are of eq. width, this is the width of every group in bp
        self.linkCount = len(self.owner.links)
        self.textMSG = txt
        self.positiveGroups = DinamicList([], 'LengthGroups:positiveGroups:')
        self.badLinkLengthCount = 0
        #self.initGroupsEqualWidth()
        self.initGroups()
        iii=9
    
    def initGroupsEqualWidth(self):
        #Length groups have the same size in bp
        lenghtMap = {} #garumam piekārto linku indeksu sarakstu ar tieši tādu (normētu) garumu
        for i in range(self.linkCount):
            v = self.owner.normLinkLength(i)
            if v not in lenghtMap:
                lenghtMap[v] = []
            lenghtMap[v].append(i)
        g = 0
        self.groupElem = {0: []} #grupas linku saraksts
        self.elemInd = [0 for i in range(self.linkCount)] #Kurā vietā links atrodas iekš groupElem
        self.elemGroup = [-1 for _ in range(self.linkCount)] #linka indeksam piekārto grupas numuru
        self.minNormLengthInGroup = [0] #katrai grupai piekārto mazāko normalizēto garumu, kas pieder šai grupai
        groupLengthLim = self.groupWidth
        for normalizedLength in sorted(lenghtMap.keys()):
            if normalizedLength>groupLengthLim:
                g+=1
                self.groupElem[g] = []
                self.minNormLengthInGroup.append(groupLengthLim+1)
                groupLengthLim += self.groupWidth
            
            for linkInd in lenghtMap[normalizedLength]:
                self.elemInd[linkInd] = len(self.groupElem[g])
                self.groupElem[g].append(linkInd)
                self.elemGroup[linkInd] = g

        self.maxExistingLength = normalizedLength #longest (normalized) link length in a given graph
        #self.minNormLengthInGroup.pop()
        self.minNormLengthInGroup.append(sys.maxsize)
        self.groupElem[g+1] = []
        self.firstSwappedLink = {g: len(self.groupElem[g]) for g in self.groupElem.keys()} # Mēs nolēmām groupElem sarakstu sākumā glabāt nekad neswapotos linku indeksus, un beigās - swapotos; 
                                                                                            #šis ir norādes uz pirmo swapoto linku šajos sarakstos
        self.groupInitLength = {g: len(self.groupElem[g]) for g in self.groupElem.keys()}
        self.groupCount = len(self.groupElem)

        self.firstGroupArea = self.groupInitLength[0]*self.minNormLengthInGroup[1] #Normazlization coeff. - area of the first, hopefully biggest group
        self.groupKeys = list(self.groupElem.keys())

    def initGroups(self):
        lenghtMap = {} #garumam piekārto linku indeksu sarakstu ar tieši tādu (normētu) garumu
        for i in range(self.linkCount):
            v = self.owner.normLinkLength(i)
            if v not in lenghtMap:
                lenghtMap[v] = []
            lenghtMap[v].append(i)
        ecount = 0
        lcount = 0
        g = 0
        self.groupElem = {0: []} #grupas linku saraksts
        self.elemInd = [0 for i in range(self.linkCount)] #Kurā vietā links atrodas iekš groupElem
        self.elemGroup = [-1 for _ in range(self.linkCount)] #linka indeksam piekārto grupas numuru
        self.minNormLengthInGroup = [0] #katrai grupai piekārto mazāko normalizēto garumu, kas pieder šai grupai
        for normalizedLength in sorted(lenghtMap.keys()):
            lcount+=1 #unique lengths in group 
            for linkInd in lenghtMap[normalizedLength]:
                self.elemInd[linkInd] = len(self.groupElem[g])
                self.groupElem[g].append(linkInd)
                ecount+=1 #edge count in group
                self.elemGroup[linkInd] = g
            if ecount >= self.minElementCount and lcount >= self.minLengthCount:
                ecount = 0
                lcount = 0
                g += 1
                self.groupElem[g] = []
                self.minNormLengthInGroup.append(normalizedLength+1)
        self.maxExistingLength = normalizedLength #longest (normalized) link length in a given graph
        #self.minNormLengthInGroup.pop()
        self.minNormLengthInGroup.append(sys.maxsize)
        self.groupElem[g+1] = []
        self.firstSwappedLink = {g: len(self.groupElem[g]) for g in self.groupElem.keys()} # Mēs nolēmām groupElem sarakstu sākumā glabāt nekad neswapotos linku indeksus, un beigās - swapotos; 
                                                                                            #šis ir norādes uz pirmo swapoto linku šajos sarakstos
        self.groupInitLength = {g: len(self.groupElem[g]) for g in self.groupElem.keys()}
        self.groupCount = len(self.groupElem)

        self.firstGroupArea = self.groupInitLength[0]*self.minNormLengthInGroup[1] #Normazlization coeff. - area of the first, hopefully biggest group
        self.groupKeys = list(self.groupElem.keys())

    def getGroupDifference(self, groupInd):
        #Knowing a group index, return number of extra links in this group, compared to number of links at the beginning
        return len(self.groupElem[groupInd])-self.groupInitLength[groupInd]

    def remove(self, linkInd):
        #linkInd - linka indekss owner.links sarakstā.
        #Pilnībā izdzēš datus par šo linku, lai pēc tam ar add to pieliktu atpakaļ
        group = self.elemGroup[linkInd] #šai garuma grupai pieder linkInd
        curIndexInGroup = self.elemInd[linkInd] #zaļā bultiņa, indekss grupas sarakstā
        lastUnswappedIndex = self.firstSwappedLink[group]-1 #pēdējais nemainītais elements sarakstā
        lastUnswappedLinkInd = self.groupElem[group][lastUnswappedIndex] #vērtība linkInd, kura bija pēdējā neswapotā
        lastLinkInd = self.groupElem[group][-1]

        if lastUnswappedIndex+1 == len(self.groupElem[group]):
            #vēl neviens nekad netika swapots
            if curIndexInGroup==len(self.groupElem[group])-1:
                #removing the last element
                self.firstSwappedLink[group]-=1
            else:
                #removing an element in between
                self.groupElem[group][curIndexInGroup] = self.groupElem[group][-1]
                self.elemInd[self.groupElem[group][-1]] = curIndexInGroup
                self.firstSwappedLink[group]-=1
        else:
            #Grupā ir vismaz 1 elements, kas jau ir ticis swapots
            if curIndexInGroup > lastUnswappedIndex:
                #mat ārā jau kādreiz noswapoto linku
                if curIndexInGroup==len(self.groupElem[group])-1:
                    #ja ir pēdējais elements
                    pass
                else:
                    self.groupElem[group][curIndexInGroup] = self.groupElem[group][-1]
                    self.elemInd[self.groupElem[group][-1]] = curIndexInGroup
            else:
                #met ārā vēl nekad neswapoto linku
                if curIndexInGroup==lastUnswappedIndex:
                    #met ārā pēdējo nekad neswapoto
                    self.groupElem[group][curIndexInGroup] = self.groupElem[group][-1]
                    self.elemInd[self.groupElem[group][-1]] = curIndexInGroup
                    self.firstSwappedLink[group]-=1
                else:
                    self.groupElem[group][curIndexInGroup] = self.groupElem[group][lastUnswappedIndex]
                    self.elemInd[self.groupElem[group][lastUnswappedIndex]] = curIndexInGroup

                    self.groupElem[group][lastUnswappedIndex] = self.groupElem[group][-1]
                    self.elemInd[self.groupElem[group][-1]] = lastUnswappedIndex

                    self.firstSwappedLink[group]-=1

        groupDiff = self.getGroupDifference(group)
        if groupDiff <= 0:
            self.badLinkLengthCount += 1
        else:
            if groupDiff == 1:
                # print("About to remove from positive group", linkGroup, self.positiveGroups, self.groupInitLength)
                self.positiveGroups.removeValue(group)  
                # print("Removed from positive group", linkGroup, self.positiveGroups, self.groupInitLength)      
            self.badLinkLengthCount -= 1

        self.groupElem[group].pop()
        self.elemInd[linkInd] = -1
        self.elemGroup[linkInd] = -1

        


    
    def getGroupForLinkInd(self, linkInd):
        return self.getGroup(self.owner.normLinkLength(linkInd))

    def getGroup(self, length):
        #length must be already normalized
        def binSearch(L, target):
            #return index of target in list L in log(n) time, target need not be in L
            #find index in L that is closest but >= than element that is the target
            n = len(L)
            if n<3 or L[0]==target: return 0
            left = 0
            right = n-1
            while right-left>1:
                midInd = (left+right)//2
                if L[midInd]==target: 
                    return midInd
                if L[midInd]<target:
                    left=midInd
                else:
                    right=midInd  
            return left

        return binSearch(self.minNormLengthInGroup, length)     
        #knowing link index, get group linearily
        # normLinkLength = self.owner.normLinkLength(linkInd)
        # answ = 0
        # trueCount = len(self.minNormLengthInGroup)
        # while answ<trueCount and (self.minNormLengthInGroup[answ])<=normLinkLength:
        #     answ+=1
        # #if answ==trueCount: return answ
        # return answ-1

    def add(self, linkInd):
        #add a link to one of the lengthGroups
        group = self.getGroupForLinkInd(linkInd)
        groupDiff = self.getGroupDifference(group)
        if groupDiff < 0:
            self.badLinkLengthCount -= 1
        else:
            if groupDiff==0:
                self.positiveGroups.addValue(group)  
            self.badLinkLengthCount += 1
        self.elemInd[linkInd] = len(self.groupElem[group])
        self.groupElem[group].append(linkInd)
        self.elemGroup[linkInd] = group

    def __getitem__(self, groupIndex):
        return self.groupElem[groupIndex]
    
    def positiveCandidate(self, chooseUnswapped=False):
        positiveGroups = self.positiveGroups.data
        if len(positiveGroups)==0:
            return self.getRandomElemFromRandomAnyGroup()
        groupWeights = [len(self.groupElem[pg]) for pg in positiveGroups]
        candGroup = random.choices(positiveGroups, weights=groupWeights, k=1)[0]
        #candGroup = self.positiveGroups.data[random.randrange(0, len(self.positiveGroups.data))]
        rangeCount = len(self.groupElem[candGroup])
        if self.firstSwappedLink[candGroup]>0:
            rangeCount = self.firstSwappedLink[candGroup]
        candInd = self.groupElem[candGroup][random.randrange(0, rangeCount)]
        #candInd = self.groupElem[candGroup][random.randrange(0, len(self.groupElem[candGroup]))]
        return candInd
    
    def getRandomElemFromRandomAnyGroup(self):
        #ņem no nejaušas grupas nejaušu linku
        return random.randrange(0, self.owner.linkCount)

    def getPositiveElemFromRandomPositiveGroup(self):
        #ņem nemainītu linku no pozitīvas grupas nejaušu linku
        positiveGroups = [pg for pg in self.positiveGroups.data if self.firstSwappedLink[pg]>0]
        return self.getRandomElemFromWeightedGroups(positiveGroups)
        
    def getRandomElemFromRandomWeightedGroup(self):
        #choose just any unchanged link, respecting weights - biežāk ņemot linkus no grupām, kuros ir relatīvi maz samainīto linku
        groups = [pg for pg in self.groupKeys if self.firstSwappedLink[pg]>0]
        return self.getRandomElemFromWeightedGroups(groups)
    
    def getRandomElemFromWeightedGroups(self, groups):
        #generic function that takes groups as agument
        if len(groups)==0: return self.getRandomElemFromRandomAnyGroup()
        groupWeights = [self.firstSwappedLink[pg]/(self.groupInitLength[pg]) for pg in groups]
        candGroup = random.choices(groups, weights=groupWeights, k=1)[0]         
        chosenLinkInd = self.groupElem[candGroup][random.randrange(0, self.firstSwappedLink[candGroup])]
        return chosenLinkInd

    
    def countBadLinks(self):
        #runs through all links and manually counts bad link count
        badLinkCount = 0
        for groupInd in self.groupElem.keys():
            badLinkCount+=abs(self.getGroupDifference(groupInd))
        return badLinkCount
    
    def getAllLinksFromPositiveGroups(self):
        #izskrien cauri visām pozitīvām grupām, uztaisa sarakstu ar visiem indeksiem ar linkiem ar pozitīviem garumiem;
        #Atgriež sarakstu ar indeksiem
        positiveGroupInds = self.positiveGroups.data #list of groupInds
        linksFromPositiveGroups = []
        for groupInd in positiveGroupInds:
            linksFromPositiveGroups.extend(self.groupElem[groupInd])
        return linksFromPositiveGroups







            







        




# class LengthGroups:
#     # Each element of data supposed to be uniq.
#     def __init__(self, owner, linkCount, minElementCount=16, minLengthCount=4, txt='Unknown'):
#         self.owner = owner
#         self.minElementCount = minElementCount
#         self.minLengthCount = minLengthCount
#         self.textMSG = txt
#         self.linkCount = linkCount
#         self.positiveGroups = DinamicList([], 'LengthGroups:positiveGroups:')
#         self.badLinkLengthCount = 0
#         self.initGroups()
#     def initGroups(self):        
#         lenghtMap = {} #garumam piekārto linku sarakstu ar tieši tādu garumu
        
#         for i in range(self.linkCount):
#             v = self.owner.normLinkLength(i)
#             if v not in lenghtMap:
#                 lenghtMap[v] = []
#             lenghtMap[v].append(i)
#         ecount = 0
#         lcount = 0
#         g = 0
#         k = 0
#         self.elemInd = [0 for i in range(self.linkCount)]
#         #self.lenGroup
#         self.lenGroup = [-1 for i in range(self.owner.maxModInd)]
#         self.lenGroup[0] = g
#         self.groupElem = {0: []}

#         for i in range(self.owner.maxModInd): #i is a possible length
#             self.lenGroup[i] = g #Garumam piekārto grupu
#             if i not in lenghtMap: #there is not a single link with this length
#                 self.lenGroup[i] = g
#             else:
#                 lcount+=1
#                 for linkInd in lenghtMap[i]:
#                     self.elemInd[linkInd] = len(self.groupElem[g])
#                     self.groupElem[g].append(linkInd)
#                     ecount+=1
#                 if ecount >= self.minElementCount and lcount >= self.minLengthCount:
#                     ecount = 0
#                     lcount = 0
#                     g += 1
#                     self.groupElem[g] = []
#         maxGroupInd = max(list(self.groupElem.keys()))
#         for g in range(maxGroupInd):
#             if g not in self.groupElem:
#                 self.groupElem[g] = []
#                 print('\ninitGroups:: create empty group', g)
#         self.groupInitLength = {g: len(self.groupElem[g]) for g in self.groupElem.keys()}
#         iii=9

#     def rollbackGroups(self, originalLengthGroup):
#         #Reset instance to another, old version of it.
#         # self.originalLenGroup = {"lenGroup": copy.deepcopy(self.lengthGroups.lenGroup), 
#         #                          "groupInitLength": copy.deepcopy(self.lengthGroups.groupInitLength)}



#         self.lenGroup = originalLengthGroup["lenGroup"]
#         self.groupInitLength = originalLengthGroup["groupInitLength"]
#         self.positiveGroups = DinamicList([], 'LengthGroups:positiveGroups:Rollback:')
#         self.groupElem = {}

#         for linkInd in range(len(self.links)): #or self.owner.links? !!!!!
#             v = self.owner.normLinkLength(linkInd)
#             vGroup = self.lenGroup[v]
#             if vGroup not in self.groupElem:
#                 self.groupElem[vGroup] = []
#             self.elemInd[linkInd] = len(self.groupElem[vGroup])
#             self.groupElem[vGroup].append(linkInd)
        
#         #Find groups that are now positive
#         for groupInd in self.groupElem.keys():
#             if len(self.groupElem[groupInd])>self.groupInitLength[groupInd]:
#                 self.positiveGroups.addValue(groupInd)

#         iii=9

#     def __str__(self):
#         S = ""
#         for i in range(len(self.lenGroup)):
#             if i in self.groupElem:
#                 v = [[j, self.linkLength(j)] for j in self.groupElem[i]]
#                 S+="\n"+ str(self.badLinkLengthCount) +str(i)+str(v)
#         S = S+"\n\n\n"
#         return S

#     def linkLength(self, linkIndex): # return length of the links
#         return self.owner.normLinkLength(linkIndex)
#     def linkGroupLength(self, linkIndex): # return length of the links length group
#         return len(self.lenGroup[self[linkIndex]])
#     def linkGroupLengthDiff(self, group): # return diff of length of current group length and the initial 
#         # len(self.groupElem[linkGroup]),
#         #           self.groupInitLength[linkGroup]
#         if group not in self.groupElem or group not in self.groupInitLength:
#             if group not in self.groupElem and group not in self.groupInitLength:
#                 return 0
#             elif group not in self.groupElem:
#                 return -self.groupInitLength[group]
#             else:
#                 return len(self.groupElem[group])
#         return len(self.groupElem[group]) - self.groupInitLength[group]
#     def __getitem__(self, linkIndex): #knowing linkIndex, get the group
#         return self.lenGroup[self.linkLength(linkIndex)]
    
#     def swapLinksInGroupWithLast(self, linkIndex):
#         candGroup = self[linkIndex] #the group of this link
#         lastInCandGroup = self.groupElem[candGroup][-1]
#         self.groupElem[candGroup][self.elemInd[linkIndex]],self.groupElem[candGroup][-1] = \
#             self.groupElem[candGroup][-1],self.groupElem[candGroup][self.elemInd[linkIndex]]
#         self.elemInd[linkIndex], self.elemInd[lastInCandGroup] = self.elemInd[lastInCandGroup], self.elemInd[linkIndex]

#     def anyCandidate(self, linkIndex):
#         candGroup = self[linkIndex] #the group of this link
#         self.swapLinksInGroupWithLast(linkIndex)

#         if len(self.groupElem[candGroup])<2:
#             return -1
#         candInd = self.groupElem[candGroup][random.randrange(0, len(self.groupElem[candGroup])-1)]
#         return candInd
    
    # def positiveCandidate(self, linkIndex):
    #     if linkIndex==-1:
    #         # if len(self.positiveGroups.data)==0:
    #         #     return
    #         candGroup = self.positiveGroups.data[random.randrange(0, len(self.positiveGroups.data))]
    #         candInd = self.groupElem[candGroup][random.randrange(0, len(self.groupElem[candGroup]))]
    #         return candInd
    #     else:
    #         return self.anyCandidate(linkIndex)

    
#     def remove(self, linkIndex):
#         #knowing a link index, remove it from groupElem, update the elemInd
#         linkGroup = self[linkIndex]

#         if linkIndex != self.groupElem[linkGroup][-1]:
#             self.groupElem[linkGroup][self.elemInd[linkIndex]] = self.groupElem[linkGroup][-1]
#             self.elemInd[self.groupElem[linkGroup][-1]] = self.elemInd[linkIndex]
#         self.elemInd[linkIndex] = -1
#         groupDiff = self.linkGroupLengthDiff(linkGroup)
#         if groupDiff <= 0:
#             self.badLinkLengthCount += 1
#         else:
#             if groupDiff == 1:
#                 # print("About to remove from positive group", linkGroup, self.positiveGroups, self.groupInitLength)
#                 self.positiveGroups.removeValue(linkGroup)  
#                 # print("Removed from positive group", linkGroup, self.positiveGroups, self.groupInitLength)      
#             self.badLinkLengthCount -= 1
        
#         self.groupElem[linkGroup].pop()
#         # print('rem::', linkIndex, linkGroup, self.badLinkLengthCount, self.groupElem[linkGroup])
    
#     def add(self, linkIndex):
#         linkGroup = self[linkIndex]
#         groupDiff = self.linkGroupLengthDiff(linkGroup)
#         if groupDiff < 0:
#             self.badLinkLengthCount -= 1
#         else:
#             if self.linkGroupLengthDiff(linkGroup) == 0:
#                 # print("About to add to positive group", linkGroup, self.positiveGroups, self.groupInitLength)  
#                 self.positiveGroups.addValue(linkGroup)  
#                 # print("Added to positive group", linkGroup, self.positiveGroups, self.groupInitLength)            
#             self.badLinkLengthCount += 1
#         self.elemInd[linkIndex] = len(self.groupElem[linkGroup])
#         self.groupElem[linkGroup].append(linkIndex)
        


















class SegmentGroups:
    # Each element of data supposed to be uniq.
    def __init__(self, owner, minElementCount=512, minSegmentCount=32, txt='Unknown'):
        self.minElementCount = minElementCount
        self.minSegmentCount = minSegmentCount #minLengthCount --> minSegmentCount
        self.textMSG = txt
        self.owner = owner
        self.positiveGroups = DinamicList([], 'SegmentGroups:')
        self.initGroups(self.owner.segments)
    def createSegmMap(self):        
        segmentMap = {}
        for link in self.owner.links:
            for i in [0, 1]:
                v = (link[i] << 1) + i
                if v not in segmentMap: segmentMap[v] = 0
                segmentMap[v] += 1
        return segmentMap
    def initGroups(self, data):        
        segmentMap = self.createSegmMap()

        self.elemInd = {}
        # maxSegment= self.owner.segments[-1]
        ecount = 0
        lcount = 0
        group = 0
        k = 0
        self.elemInd = [-1 for i in range(len(self.owner.links) * 2 + 1)]
        self.segGroup = {}
        for v in self.owner.segments:
            self.segGroup[v << 1] = -1
            self.segGroup[(v << 1) + 1] = -1
        self.groupElem = {0: []}
        segmentKeys = sorted(list(self.segGroup.keys()))
        for i in range(len(segmentKeys)):
            self.segGroup[i] = group
            if i in segmentMap.keys():
                lcount+=1
                if i not in segmentMap:
                    continue
                ecount += segmentMap[i]
                if ecount >= self.minElementCount and lcount >= self.minSegmentCount:
                    ecount = 0
                    lcount = 0
                    group += 1
                    # self.groupElem[group] = []
        iii=9
        del segmentMap
        for linkIndex in range(len(self.owner.links)):
            self.add(linkIndex)
        self.groupElemMaxKey = max(self.groupElem.keys())+1

    def __getitem__(self, index): #knowing linkIndex, get the group
        return self.segGroup[self.owner.links[index >> 1][index & 1]]

    def remove(self, linkIndex):
        #knowing a link index, remove it from groupElem, update the elemInd
        for i in [0,1]:
            index = (linkIndex << 1) + i
            group = self[index]
            if index != self.groupElem[group][-1]:
                self.groupElem[group][self.elemInd[index]] = self.groupElem[group][-1]
                self.elemInd[self.groupElem[group][-1]] = self.elemInd[index]
            self.elemInd[index] = -1
            self.groupElem[group].pop()
        iii=9
    
    def add(self, linkIndex):
        for i in [0,1]:
            index = (linkIndex << 1) + i
            group = self[index]
            if group not in self.groupElem:
                self.groupElem[group] = []
            self.elemInd[index] = len(self.groupElem[group])
            self.groupElem[group].append(index)
    
    def isValid(self, linkIndex):
        for i in [0, 1]:
            index = (self.owner.links[linkIndex][i] << 1) + i
            group = self[index]
            if group not in self.segGroup:
                warningMSG(self.textMSG + '::isValid:: link has not group', [i, self.owner.links[linkIndex][:2], index])
            if self.groupElem[group][self.elemInd[index]] != index:
                warningMSG(self.textMSG + '::isValid:: link endpoint is not in group', [i, self.owner.links[linkIndex][:2], index])
        return True
    
    def isAllValid(self):
        for linkIndex in range(len(self.owner.links)):
            for i in [0, 1]:
                index = (self.owner.links[linkIndex][i] << 1) + i
                group = self[index]
                if group not in self.segGroup:
                    warningMSG(self.textMSG + '::isAllValid:: link has not group', [i, self.owner.links[linkIndex][:2], index, group])
                if self.groupElem[group][self.elemInd[index]] != index:
                    warningMSG(self.textMSG + '::isAllValid:: link endpoint is not in group', [i, self.owner.links[linkIndex][:2], index])
        indexSet = None
        icount = 0
        for group in self.groupElem.keys():
            gset = set(self.groupElem[group])
            icount += len(self.groupElem[group])
            if len(gset) != len(self.groupElem[group]): 
                warningMSG(self.textMSG + '::isAllValid:: duplicate index in self.segGroup[group]', [group, self.links[linkIndex][:2], index, self.groupElem[group]])
            if indexSet is None:
                indexSet = gset
            elif len(gset) > 0:
                indexSet = indexSet.union(gset)
        if len(indexSet) != len(self.owner.links) * 2:
            warningMSG(self.textMSG + '::isAllValid:: len(indexSet) != len(self.links) * 2', [len(indexSet), len(self.owner.links) * 2])
        return True

