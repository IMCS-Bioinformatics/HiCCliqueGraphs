#!/usr/bin/env python
import csv, os, math, time, copy
import re
# import Utils
import CReader

def perc(a, b):
    if b == 0: return 0
    return round(a / b, 6)

class Graph():
    def __init__(self, ch, data, minBase): 
        self.ch = ch
        self.minBase = minBase
        links = data['chrValues'][ch]['links']
        self.links = [links[i][:3] for i in range(len(links))]
        self.tissueBits = data['tissueBits']
        self.tissues = data['tissueIDs']
        self.n = len(data['chrValues'][ch]['segments'])
        self.ckeys = ['links', 'triangles', 'bases']        
        self.rez = {}
        self.sumValues = {}
        self.sumCount = 0
        m = len(self.tissues)
        for i in range(m - 1):
            for j in range(i, m):
                self.sumValues[(self.tissues[i], self.tissues[j])] = [[0, 0, -1] for k in self.ckeys]
        self.start_time = time.time()
        self.initData()
        iii = 1
    def getTime(self):
        ctime = time.time() - self.start_time
        self.start_time = time.time()
        return ctime
    def setAdj(self, links):
        adj = [set() for i in range(self.n)]
        for e in self.links:
            adj[e[0]].add(e[1])
        return adj
    def initData(self):
        # self.adj = [set() for i in range(self.n)]
        self.inDeg = [0 for i in range(self.n)]
        self.outDeg = [0 for i in range(self.n)]
        self.linkTissues = {}
        for e in self.links:
            # self.adj[e[0]].add(e[1])
            self.outDeg[e[0]] += 1
            self.inDeg[e[1]] += 1
            self.linkTissues[(e[0], e[1])] = e[2]
        self.segmDeg = [self.inDeg[i] + self.outDeg[i] for i in range(self.n)]
        self.segmOvelap = [0 for i in range(self.n)]
        s = 0
        for i in range(self.n):
            self.segmOvelap[i] = s
            s += self.outDeg[i] - self.inDeg[i]
        self.segmOvelap[self.n - 1] = s
        self.calcTriangles()
        self.calcBases()
        iii = 1
    def calcTriangles(self):
        self.triangles = []
        adj = self.setAdj(self.links)
        for e in self.links:
            comm = adj[e[0]].intersection(adj[e[1]])
            triangles = []
            for v in comm:
                bt = e[2] & self.linkTissues[(e[0], v)] & self.linkTissues[(e[1], v)]
                if bt > 0:
                    triangles.append((e[0], e[1], bt))
            self.triangles += triangles
        self.triangles.sort(key=lambda x: (x[0], x[1], x[2]))
        del adj
    def calcBases(self):
        bases = {}
        for t in self.tissueBits.keys():
            bt = self.tissueBits[t]
            links = [v for v in self.links if bt & v[2] == bt]
            adj = self.setAdj(links)
            for e in links:
                comm = adj[e[0]].intersection(adj[e[1]])
                if len(comm) >= self.minBase:
                    v = (e[0], e[1])
                    if v not in bases: bases[v] = 0
                    bases[v] |= bt
        self.bases = [[v[0], v[1], bases[v]] for v in bases.keys()]
        self.bases.sort(key=lambda x: (x[0], x[1], x[2]))
    def outData(self):
        return { \
            'inDeg': self.inDeg, \
            'outDeg': self.outDeg, \
            'segmDeg': self.segmDeg, \
            'segmOvelap': self.segmOvelap, \
            'segmOvelap': self.segmOvelap, \
            'triangles': self.triangles, \
            'bases': self.bases, \
            'links': self.links, \
            'tissueBits': self.tissueBits, \
            'ch': self.ch,\
        }

class Compare():
    def __init__(self, bdata, tdata, tisBits): 
        header = ['chr', 'tissue1', 'tissue2', 'comm links', 'test links', 'orig links', 'comm tri', 'test tri', 'orig tri', 'comm base', 'test base', 'orig base']
        self.csvRez = [header]
        #for ch in bdata.keys():
        for ch in [bdata["ch"]]:
            print('\t', ch)
            self.baseData = bdata
            self.testData = tdata
            self.tisBits = tisBits
            self.tissues = list(self.tisBits.keys())
            self.ckeys = ['links', 'triangles', 'bases']
            m = len(self.tissues)
            for i in range(m - 1):
                for j in range(i, m):
                    values = [ch, self.tissues[i], self.tissues[j]]
                    for k in range(len(self.ckeys)):
                        key = self.ckeys[k]
                        comm = self.testCommon(key, self.tissues[i], self.tissues[j])
                        values = values + comm[:3]
                    self.csvRez.append(values)
                iii = 1 
            iii = 1 
        iii = 1
    def getTime(self):
        ctime = time.time() - self.start_time
        self.start_time = time.time()
        return ctime
    def testCommon(self, key, t0, t1):
        bt = self.tisBits[t0] + self.tisBits[t1]
        baseSet = set([tuple(v) for v in self.baseData[key] if (v[2] & bt) == bt])
        testSet = set([tuple(v) for v in self.testData[key] if (v[2] & bt) == bt])
        comData = list(baseSet.intersection(testSet))
        comData.sort(key=lambda k: (k[0], k[1], k[2]))
        # print('\t\t\t', key, len(comData), len(baseSet), len(testSet), perc(len(comData), len(baseSet)))
        commonData = [list(v) for v in comData]
        return [len(comData), len(testSet), len(baseSet), commonData]


import json
fn = "sampleData/data-pvalue-5-fin-min.json" #<<<< Original graph
f = open(fn)
dataOrig = json.load(f) #list of links and segments
f.close()


fn = "sampleData/data-pvalue-5-fin-minRND.json"
f = open(fn)
dataRand = json.load(f) #list of links and segments
f.close()






Grand = Graph(ch="chr6", data=dataRand, minBase=16).outData()
Gorig = Graph(ch="chr6", data=dataOrig, minBase=16).outData()
C = Compare(bdata=Gorig, tdata=Grand, tisBits=dataOrig["tissueBits"])

import csv
with open('RealVsRandCompareResults.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerows(C.csvRez)

    for c in [15]: #chromosomes to compare
        ch = "chr"+str(c)
        Grand = Graph(ch=ch, data=dataRand, minBase=16).outData()
        Gorig = Graph(ch=ch, data=dataOrig, minBase=16).outData()
        C = Compare(bdata=Gorig, tdata=Grand, tisBits=dataOrig["tissueBits"])
        spamwriter.writerows(C.csvRez[1:])
    
    # ch = "chrX"
    # Grand = Graph(ch=ch, data=dataRand, minBase=16).outData()
    # Gorig = Graph(ch=ch, data=dataOrig, minBase=16).outData()
    # C = Compare(bdata=Gorig, tdata=Grand, tisBits=dataOrig["tissueBits"])
    # spamwriter.writerows(C.csvRez[1:])