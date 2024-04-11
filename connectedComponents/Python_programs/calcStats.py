#!/usr/bin/python
import math
import CReader
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--pvalue', type=float, default=6)
parser.add_argument('--dataPath', default=os.path.join(os.path.dirname(
    __file__), "../Processed_data/Normal_Hi-C/"))
args = parser.parse_args()

pvalue = args.pvalue
if not isinstance(pvalue, int) and pvalue.is_integer():
    pvalue = int(pvalue)
inputJSON = args.dataPath + "Interactions_with_states/Pvalue" + str(
    pvalue) + "/"
inputCSV = args.dataPath + "Components/Pvalue" + str(pvalue) + "/"
inputTADs = args.dataPath + "TADs/"

outputCSV = args.dataPath + "Components_extra/Pvalue" + str(pvalue) + "/"


def connComp(v, w, links):
    ccomp = set()
    adj = {}
    for e in links:
        if e[0] not in adj:
            adj[e[0]] = []
        adj[e[0]].append(e[1])
        if e[1] not in adj:
            adj[e[1]] = []
        adj[e[1]].append(e[0])
    stack = [v]
    ccomp.add(v)
    while len(stack) > 0:
        v = stack.pop()
        for w in adj[v]:
            if w not in ccomp:
                stack.append(w)
                ccomp.add(w)
    return sorted(ccomp)


def hasProtein(a, geneSymbolsProteins):
    for gene in a["genes"]:
        if geneSymbolsProteins[gene]["protein"] != "":
            return True
    return False


def getBitCount(bits):
    i = 1
    count = 0
    while i <= bits:
        if (bits & i) == i:
            count += 1
        i *= 2
    return count


chrs = [str(i) for i in range(1, 23)] + ["X"]
chrs.reverse()
for chr in chrs:
    data = CReader.readJson(inputJSON + "chr" + chr + ".json")
    tads = 0
    if os.path.exists(inputTADs + "RealRes-chr" + chr + "-SpectralTad.json"):
        tads = CReader.readJson(
            inputTADs + "RealRes-chr" + chr + "-SpectralTad.json")
    components = CReader.readCsvToDic(inputCSV + "chr" + chr + ".csv")

    prev_header = CReader.readCsvHeader(inputCSV + "chr" + chr + ".csv")
    header = prev_header + ['geneCount',
                            'proteinCount', 'minStateCount']
    if tads:
        header = prev_header + ['geneCount',
                                'proteinCount', 'averageTads', 'minStateCount', 'maxTadGap']
    rez = [header]
    allLinks = data["links"]
    segments = [(v["segment"][0] + v["segment"][1]) //
                2 for v in data["segments"]]
    nodes = [v for v in data["segments"]]
    geneSymbolsProteins = data["geneSymbolsProteins"]
    k = 0
    for comp in components:
        k += 1
        bt = int(comp["basetissues"])
        compLinks = [[v["bin1"], v["bin2"]]
                     for v in allLinks if (v["bits"] & bt) == bt]
        ccomp = connComp(int(comp["from"]), int(comp["to"]), compLinks)
        gene = [v for v in ccomp if len(nodes[v]["genes"])]
        prot = [v for v in ccomp if hasProtein(nodes[v], geneSymbolsProteins)]
        baseTissues = [tissue for tissue in data["tissueData"]["tissueBits"] if (
            bt & data["tissueData"]["tissueBits"][tissue]) == data["tissueData"]["tissueBits"][tissue]]
        if tads:
            tadSum = 0
            maxTadGap = 0
            for tissue in tads["TADs"].keys():
                i = 0
                tadSet = set()
                for v in ccomp:
                    while i+1 < len(tads["TADs"][tissue]) and tads["TADs"][tissue][i+1] <= segments[v]:
                        i += 1
                        if i == len(tads["TADs"][tissue]):
                            print("####", tissue, segments[v], tads[tissue])
                    tadSet.add(i)
                if tissue in baseTissues:
                    tadList = sorted(list(tadSet))
                    for i in range(len(tadList)-1):
                        if tadList[i+1] - tadList[i] > 1:
                            maxTadGap = max(
                                maxTadGap, tadList[i+1] - tadList[i] - 1)
                tadSum += len(tadSet)
            tadAver = round(tadSum / len(tads["TADs"]), 2)
            comp['averageTads'] = tadAver
            comp['maxTadGap'] = maxTadGap
        minStateCount = (1 << 31)
        for tissue in baseTissues:
            if tissue in data["tissueData"]["tissuesWithStates"]:
                states = set()
                for v in ccomp:
                    if tissue in nodes[v]["states"]:
                        for state in nodes[v]["states"][tissue]:
                            states.add(state)
                minStateCount = min(len(states), minStateCount)
        # If none of the base tissues has state data
        if minStateCount == (1 << 31):
            minStateCount = -1

        comp['geneCount'] = len(gene)
        comp['proteinCount'] = len(prot)
        comp['minStateCount'] = minStateCount
        values = [comp[k] for k in header]
        rez.append(values)
    if not os.path.exists(outputCSV):
        os.makedirs(outputCSV)
    CReader.saveCsv(outputCSV + "chr" + chr + ".csv", rez)
    del data
    del tads
    del components
    del allLinks
    del segments
    del nodes
