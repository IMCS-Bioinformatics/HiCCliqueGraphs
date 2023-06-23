# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 09:58:08 2023
Recalculated some GOEA analyzes, using QD classes
@author: asizo
fns can be changed
!Note, script should be run from root directory of the repository (From the directory which has the sampleData subdirectory)
"""

from universal import *
from topologicalFeatures import *
from GOEAa import *



fns = ["sampleData\\data-pvalue-5-fin-min.json","sampleData\\data-pvalue-5-fin-minRND.json",
       "sampleData\\data-pvalue-0.7-fin-min.json","sampleData\\data-pvalue-0.7-fin-minRND.json",
       "sampleData\\data-pvalue-10-fin-min.json","sampleData\\data-pvalue-10-fin-minRND.json"
    ]

Us = [UniversalDS(fn) for fn in fns]
Us[0].DS, Us[1].DS = "Blood cell pcHi-C", "Randomized blood cell pcHi-C"
Us[2].DS, Us[3].DS = "Tissue pcHi-C", "Randomized tissue pcHi-C"
Us[4].DS, Us[5].DS = "Tissue Hi-C", "Randomized tissue Hi-C"

Biolog = GOEATool("GOEA/ontologyData/go-basic.obo", "GOEA/ontologyData/hg19_genes_w-go.txt")


#triangles1 vs links1
rezFN = "CL1vsLinks1-{}vsRand.csv".format(Us[0].DS)
what = "Cycles C3 w/ 1+ tissue vs Links w/ 1+ tissue"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=1)
        C = Cliques(chData, minC3TissueCount=1)
        Biolog.GOEA(C, chData, what=what)
        print("*************>>>> new count",len(Biolog.curAccumulatedResults))
        print(Biolog.curAccumulatedResults[0:min(len(Biolog.curAccumulatedResults),3)])

Biolog.dump(rezFN)
Biolog.reset()

#triangles2 vs links2
rezFN = "CL2vsLinks2-{}vsRand.csv".format(Us[0].DS)
what = "Cycles C3 w/ 2+ tissue vs Links w/ 2+ tissue"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=2)
        C = Cliques(chData, minC3TissueCount=2)
        Biolog.GOEA(C, chData, what=what)
        print("*************>>>> new count",len(Biolog.curAccumulatedResults))
        print(Biolog.curAccumulatedResults[0:min(len(Biolog.curAccumulatedResults),3)])

Biolog.dump(rezFN)
Biolog.reset()
##################################################################

#triangles1 vs all possibles
rezFN = "CL2vsAll-{}vsAllPossible.csv".format(Us[0].DS)
what = "Cycles C3 w/ 2+ tissue vs All possible genes"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch)
        C = Cliques(chData, minC3TissueCount=2)
        Biolog.GOEA(C, popObject="all", what=what)
        print("*************>>>> new count",len(Biolog.curAccumulatedResults))
        print(Biolog.curAccumulatedResults[0:min(len(Biolog.curAccumulatedResults),3)])

Biolog.dump(rezFN)
Biolog.reset()

##################################################################
#Triangles 5 vs links 1
rezFN = "CL5vsLinks2.csv"
what = "Cliques w/ 5+ tissue vs Links w/ 2+ tissue"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=2)
        C = Cliques(chData, minC3TissueCount=5)
        Biolog.GOEA(C, chData, what=what)

Biolog.dump(rezFN)
Biolog.reset()

##########################################################
#links 5 vs links 1
rezFN = "Links5vsLinks1.csv"
what = "Links w/ 5+ tissues vs Links w/ 1+ tissue"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=1)
        C = ChrData(owner=U, ch=ch, minLinkTissueCount=5)
        Biolog.GOEA(C, chData, what=what)

Biolog.dump(rezFN)
Biolog.reset()
#############################################################
#bases 5 vs links 2
rezFN = "BasesLogvsLinks2-{}vsRand.csv".format(Us[0].DS)
what = "bases deg16 vs Links w/ 2+ tissue"
import numpy as np
for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=2)
        C = BasesOfBases(owner=chData)
        log = int(np.round(np.log2(len(chData.links))))
        C.reduce(deg=log)
        Biolog.GOEA(C, chData, what="bases deg{log} vs Links w/ 2+ tissue".format(log=log))

Biolog.dump(rezFN)
Biolog.reset()
################################################################
# #bases log(linkCount) vs links 1
rezFN = "BasesLogvsLinks1.csv"
what = "bases degLog vs Links w/ 1+ tissue"

for U in Us:
    for ch in U.chrs:
        print("<<<<<<<>>>>>>>", U.fn, ch)
        chData = ChrData(owner=U, ch=ch, minLinkTissueCount=1)
        C = BasesOfBases(owner=chData)
        deg = int(np.log2(len(C.links)))
        C.reduce(deg=deg)
        Biolog.GOEA(C, chData, what=f"bases deg {deg} (log link count) vs Links w/ 1+ tissue")

Biolog.dump(rezFN)
Biolog.reset()
################################################
#bases 5 vs cliques 1
# rezFN = "BasesLogvsCL2.csv"
# what = "bases degLog vs Cliques w/ 1+ tissue"

# for U in Us:
#     for ch in U.chrs:
#         print("<<<<<<<>>>>>>>", U.fn, ch)
#         chData = ChrData(owner=U, ch=ch, minLinkTissueCount=1)
#         C = BasesOfBases(owner=chData)
#         log = int(np.round(np.log2(len(chData.links))))
#         C.reduce(deg=log)
#         cliques = Cliques(chData, minC3TissueCount=1)
#         Biolog.GOEA(C, cliques, what=f"bases deg{log} vs Cliques w/ 1+ tissue")

# Biolog.dump(rezFN)
# Biolog.reset()
    
# for curFn in fns:
#     U = UniversalDS(curFn)
    




# chData = ChrData(U, ch="chr6", minLinkTissueCount=1)
# C = Cliques(chData, minC3TissueCount=2, tissueMask=0) 


# Biolog = GOEATool("C:/Users/asizo/Documents/myPrograms/go-basic.obo", "C:/Users/asizo/Documents/myPrograms/hg19_genes_w-go.txt")
# print(">>>>>>", Biolog.GOEA(stuObject=C, popObject=chData, what="Cliques w/ 2 tissues vs links w/ 1 tissue"))
# chData = ChrData(U, ch="chr19", minLinkTissueCount=1)
# C = Cliques(chData, minC3TissueCount=2, tissueMask=0) 
# Biolog.GOEA(stuObject=C, popObject=chData, what="Cliques w/ 2 tissues vs links w/ 1 tissue")
# Biolog.dump("testnr20jan.csv")