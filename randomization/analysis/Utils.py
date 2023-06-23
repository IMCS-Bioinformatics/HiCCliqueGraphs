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
    return [v for v in os.listdir(path) if v.find(ext) > 0]
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
    # Each element of data supposed to be uniq.
    def __init__(self, owner, segments, minElementCount=16, minLengthCount=4, txt='Unknown'):
        self.owner = owner
        self.minElementCount = minElementCount
        self.minLengthCount = minLengthCount
        self.textMSG = txt
        self.links = segments
        self.positiveGroups = DinamicList([], 'LengthGroups:positiveGroups:')
        self.badLinkLengthCount = 0
        self.initGroups(segments)
    def initGroups(self, data):        
        lenghtMap = {} #garumam piekārto linku sarakstu ar tieši tādu garumu
        #maxSegment = max([el[i] for el in self.links for i in [0,1]])+1
        
        for i in range(len(data)):
            #v = abs(data[i][1] - data[i][0])
            v = self.owner.normLinkLength(i)
            if v not in lenghtMap:
                lenghtMap[v] = []
            lenghtMap[v].append(i)
        maxLength = max(list(lenghtMap.keys()))*2
        ecount = 0
        lcount = 0
        g = 0
        k = 0
        # self.groupFirst = [0]
        self.elemInd = [0 for i in range(len(data))]
        self.lenGroup = [-1 for i in range(maxLength)]
        self.lenGroup[0] = g
        self.groupElem = {0: []}

        for i in range(maxLength): #i is a possible length
            self.lenGroup[i] = g
            if i not in lenghtMap: #there is not a single link with this length
                self.lenGroup[i] = g
            else:
                lcount+=1
                for linkInd in lenghtMap[i]:
                    self.elemInd[linkInd] = len(self.groupElem[g])
                    self.groupElem[g].append(linkInd)
                    ecount+=1
                if ecount >= self.minElementCount and lcount >= self.minLengthCount:
                    ecount = 0
                    lcount = 0
                    g += 1
                    # self.groupFirst.append(v) 
                    self.groupElem[g] = []
        maxGroupInd = max(list(self.groupElem.keys()))
        for g in range(maxGroupInd):
            if g not in self.groupElem:
                self.groupElem[g] = []
                print('\ninitGroups:: create empty group', g)
        self.groupInitLength = {g: len(self.groupElem[g]) for g in self.groupElem.keys()}

    # def isNegativeGroup(self, g):
    #     return self.groupInitLength[g] > len(self.groupElem[g])
    # def isNeitralGroup(self, g):
    #     return self.groupInitLength[g] == len(self.groupElem[g])

    def __str__(self):
        S = ""
        for i in range(len(self.lenGroup)):
            if i in self.groupElem:
                #v = [[j, self.links[j][1] - self.links[j][0]] for j in self.groupElem[i]]
                v = [[j, self.linkLength(j)] for j in self.groupElem[i]]
                S+="\n"+ str(self.badLinkLengthCount) +str(i)+str(v)
                #print('\n', i, v)
        S = S+"\n\n\n"
        return S

    def linkLength(self, linkIndex): # return length of the links
        #return abs(self.links[linkIndex][1]-self.links[linkIndex][0])
        return self.owner.normLinkLength(linkIndex)
    def linkGroupLength(self, linkIndex): # return length of the links length group
        return len(self.lenGroup[self[linkIndex]])
    def linkGroupLengthDiff(self, group): # return diff of length of current group length and the initial 
        # len(self.groupElem[linkGroup]),
        #           self.groupInitLength[linkGroup]
        return len(self.groupElem[group]) - self.groupInitLength[group]
    def __getitem__(self, linkIndex): #knowing linkIndex, get the group
        return self.lenGroup[self.linkLength(linkIndex)]
    
    def swapLinksInGroupWithLast(self, linkIndex):
        candGroup = self[linkIndex] #the group of this link
        lastInCandGroup = self.groupElem[candGroup][-1]
        self.groupElem[candGroup][self.elemInd[linkIndex]],self.groupElem[candGroup][-1] = \
            self.groupElem[candGroup][-1],self.groupElem[candGroup][self.elemInd[linkIndex]]
        self.elemInd[linkIndex], self.elemInd[lastInCandGroup] = self.elemInd[lastInCandGroup], self.elemInd[linkIndex]

    def anyCandidate(self, linkIndex):
        #candGroup = random.randrange(0, len(self.groupElem))
        candGroup = self[linkIndex] #the group of this link
        self.swapLinksInGroupWithLast(linkIndex)

        if len(self.groupElem[candGroup])<2:
            return -1
        candInd = self.groupElem[candGroup][random.randrange(0, len(self.groupElem[candGroup])-1)]
        return candInd
    
    def positiveCandidate(self, linkIndex):
        if linkIndex==-1:
            candGroup = self.positiveGroups.data[random.randrange(0, len(self.positiveGroups.data))]
            candInd = self.groupElem[candGroup][random.randrange(0, len(self.groupElem[candGroup]))]
            return candInd
        else:
            return self.anyCandidate(linkIndex)

    
    def remove(self, linkIndex):
        #knowing a link index, remove it from groupElem, update the elemInd
        linkGroup = self[linkIndex]
        # if linkGroup==232:
        #     iii=9
        #     print(">>>>>>>>>>rem beg", linkIndex, linkGroup, (linkGroup in self.positiveGroups.dataMap), len(self.groupElem[linkGroup]),
        #           self.groupInitLength[linkGroup], self.positiveGroups.data, self.linkGroupLengthDiff(linkGroup),
        #           self.groupElem[linkGroup])

        if linkIndex != self.groupElem[linkGroup][-1]:
            self.groupElem[linkGroup][self.elemInd[linkIndex]] = self.groupElem[linkGroup][-1]
            self.elemInd[self.groupElem[linkGroup][-1]] = self.elemInd[linkIndex]
        self.elemInd[linkIndex] = -1
        groupDiff = self.linkGroupLengthDiff(linkGroup)
        if groupDiff <= 0:
            self.badLinkLengthCount += 1
        else:
            if groupDiff == 1:
                # print("About to remove from positive group", linkGroup, self.positiveGroups, self.groupInitLength)
                self.positiveGroups.removeValue(linkGroup)  
                # print("Removed from positive group", linkGroup, self.positiveGroups, self.groupInitLength)      
            self.badLinkLengthCount -= 1
        
        # if linkGroup==232:
        #     iii=9
        #     print(">>>>>>>>>>rem end",linkIndex, linkGroup, (linkGroup in self.positiveGroups.dataMap), len(self.groupElem[linkGroup]),
        #           self.groupInitLength[linkGroup], self.positiveGroups.data, self.linkGroupLengthDiff(linkGroup),
        #           self.groupElem[linkGroup])
        self.groupElem[linkGroup].pop()
        # print('rem::', linkIndex, linkGroup, self.badLinkLengthCount, self.groupElem[linkGroup])
    
    def add(self, linkIndex):
        linkGroup = self[linkIndex]
        # if linkGroup==232:
        #     iii=9
        #     print(">>>>>>>>>>add begin ", linkIndex, linkGroup, (linkGroup in self.positiveGroups.dataMap), len(self.groupElem[linkGroup]),
        #           self.groupInitLength[linkGroup], self.positiveGroups.data, self.linkGroupLengthDiff(linkGroup), self.linkGroupLengthDiff(linkGroup),
        #           self.groupElem[linkGroup])
        
        groupDiff = self.linkGroupLengthDiff(linkGroup)
        if groupDiff < 0:
            self.badLinkLengthCount -= 1
        else:
            if self.linkGroupLengthDiff(linkGroup) == 0:
                # print("About to add to positive group", linkGroup, self.positiveGroups, self.groupInitLength)  
                self.positiveGroups.addValue(linkGroup)  
                # print("Added to positive group", linkGroup, self.positiveGroups, self.groupInitLength)            
            self.badLinkLengthCount += 1
        self.elemInd[linkIndex] = len(self.groupElem[linkGroup])
        self.groupElem[linkGroup].append(linkIndex)
        # if linkGroup==232:
        #     iii=9
        #     print(">>>>>>>>>>add end ",linkIndex, linkGroup, (linkGroup in self.positiveGroups.dataMap), len(self.groupElem[linkGroup]),
        #           self.groupInitLength[linkGroup], self.positiveGroups.data, self.linkGroupLengthDiff(linkGroup),
        #           self.groupElem[linkGroup])
        # print('add::', linkIndex, linkGroup, self.badLinkLengthCount, self.groupElem[linkGroup])

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
class SegmentGroups:
    # Each element of data supposed to be uniq.
    def __init__(self, segments, minElementCount=512, minSegmentCount=32, txt='Unknown'):
        self.minElementCount = minElementCount
        self.minSegmentCount = minSegmentCount #minLengthCount --> minSegmentCount
        self.textMSG = txt
        self.links = segments
        self.positiveGroups = DinamicList([], 'SegmentGroups:')
        self.initGroups(segments)
    def initGroups(self, data):        
        segmentMap = {}
        maxSegment= 0
        for link in self.links:
            for j in [0, 1]:
                if maxSegment < link[j]: 
                    maxSegment = link[j]
                v = (link[j]<<1) + j
                if v not in segmentMap:
                    segmentMap[v] = 0
                segmentMap[v] += 1
        ecount = 0
        lcount = 0
        group = 0
        k = 0
        segmRange = range(maxSegment*2+1)
        self.elemInd = [-1 for i in range(len(self.links) * 2 + 1)]
        self.segGroup = [-1 for i in segmRange] #index
        self.segGroup[0] = group
        self.groupElem = {0: []}

        for i in segmRange:
            self.segGroup[i] = group
            if i in segmentMap.keys():
                lcount+=1
                # self.elemInd[i] = len(self.groupElem[group])
                # self.groupElem[group].append(i)
                if i not in segmentMap:
                    continue
                ecount += segmentMap[i]
                if ecount >= self.minElementCount and lcount >= self.minSegmentCount:
                    ecount = 0
                    lcount = 0
                    group += 1
                    # self.groupElem[group] = []
        self.groupElemMaxKey = group + 1
        iii=9
        del segmentMap
        for linkIndex in range(len(self.links)):
            self.add(linkIndex)

    def __getitem__(self, index): #knowing linkIndex, get the group
        return self.segGroup[self.links[index >> 1][index & 1]]

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
            index = (self.links[linkIndex][i] << 1) + i
            group = self[index]
            if group not in self.segGroup:
                warningMSG(self.textMSG + '::isValid:: link has not group', [i, self.links[linkIndex][:2], index])
            if self.groupElem[group][self.elemInd[index]] != index:
                warningMSG(self.textMSG + '::isValid:: link endpoint is not in group', [i, self.links[linkIndex][:2], index])
        return True
    
    def isAllValid(self):
        for linkIndex in range(len(self.links)):
            for i in [0, 1]:
                index = (self.links[linkIndex][i] << 1) + i
                group = self[index]
                if group not in self.segGroup:
                    warningMSG(self.textMSG + '::isAllValid:: link has not group', [i, self.links[linkIndex][:2], index, group])
                if self.groupElem[group][self.elemInd[index]] != index:
                    warningMSG(self.textMSG + '::isAllValid:: link endpoint is not in group', [i, self.links[linkIndex][:2], index])
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
        if len(indexSet) != len(self.links) * 2:
            warningMSG(self.textMSG + '::isAllValid:: len(indexSet) != len(self.links) * 2', [len(indexSet), len(self.links) * 2])
        return True

# random.seed(12345)
# rg = [2, 1024]
# data = [sorted([random.randrange(rg[0], rg[1]), random.randrange(rg[0], rg[1])]) for i in range(1024)]
# g = SegmentGroups(data, 16, 4)
# iii = 1

# #10 reizes pamainísim linkiem galus
# n=1000
# g.isAllValid()
# for i in range(n):
#     #izvelas 2 gadijuma linkus
#     aInd = random.randrange(0, len(data))
#     bInd = random.randrange(0, len(data))
#     if aInd == bInd:
#         continue
   
#     print("Before remove", aInd, bInd)
#     g.remove(aInd)
#     g.remove(bInd)
#     # print(g)

#     data[aInd][0], data[bInd][0] = data[bInd][0], data[aInd][0]

#     g.add(aInd)
#     g.add(bInd)

#     g.isValid(aInd)
#     g.isValid(bInd)
#     g.isAllValid()

#     iii=9
# iii = 1
