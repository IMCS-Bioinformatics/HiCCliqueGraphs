import pandas as pd
import os
import json

templ = {
    "dataPath": "./sourceData/PCHiC_peak_matrix_cutoff5.txt",
    "outputPath": "./transformedData/BloodCellPCHi-C",
    "tissueIDs":  ['Mon', 'Mac0', 'Mac1', 'Mac2', 'Neu', 'MK', 'EP', 'Ery', 'FoeT', 'nCD4', 'tCD4', 'aCD4', 'naCD4', 'nCD8', 'tCD8', 'nB', 'tB'],
    "dataInds": {"firstPvalue": 11, "baitChr": 0, "baitStart": 1, "baitEnd": 2, "oeChr": 5, "oeStart": 6, "oeEnd": 7},
    "pvalue": 5
}


# Contains data for one chromosome
class DataChr:
    def __init__(self, owner):
        self.owner = owner
        self.segmCount = 0
        self.segmInd = {}
        self.linkTissues = {}
        self.segments = set()
        self.nodes = []

    # Creates an ID for segment if it does not already have one
    def indSegm(self, s):
        if s not in self.segmInd:
            self.segmInd[s] = self.segmCount
            self.segmCount += 1
        return self.segmInd[s]

    # Adds tissue to linkTissues of given link
    def setLinkInfo(self, v, tissue):
        link = tuple(sorted(v[0:2]))
        if link not in self.linkTissues:
            self.linkTissues[link] = [0, {}]
        self.linkTissues[link][0] |= self.owner.tissueBits[tissue]
        # If this link already has this tissue, choose greatest pvalue
        self.linkTissues[link][1][self.owner.tissueInd[tissue]] = max(self.linkTissues[link][1].get(self.owner.tissueInd[tissue], 0), round(
            v[3], 4))
        return link

    # Returns link as [one end ID, other end ID, tissue bitmap, tissue pvalue list]
    def getLinkRow(self, link):
        # Alert if one of ends (segments) does not have an ID
        if link[0] not in self.segmInd or link[1] not in self.segmInd:
            print("!!!Segment", link[0], "or", link[1], "without id!!!")
        u = self.indSegm(link[0])
        v = self.indSegm(link[1])
        rez = [u, v, self.linkTissues[link][0], self.linkTissues[link][1]]
        return rez


# Contains data for all chromosomes
class DataHiC:
    def __init__(self, templ):
        self.tissueIDs = sorted(templ["tissueIDs"], key=str.casefold)
        self.chrNames = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
        self.chInd = {}
        for i in range(len(self.chrNames)):
            self.chInd[self.chrNames[i]] = i
        self.inputPath = templ["dataPath"]
        self.pvalue = templ["pvalue"]
        self.tissueBits = {}
        self.tissueInd = {}
        # For each tissue assign a bit (a power of two)
        bitNames = self.tissueIDs
        for i in range(len(bitNames)):
            self.tissueBits[bitNames[i]] = 1 << i
        for i in range(len(bitNames)):
            self.tissueInd[bitNames[i]] = i
        self.chrs = [DataChr(self) for i in range(len(self.chrNames))]
        self.initData()

    # Returns needed data from a row in PCHi-C data (one end segment, other end segment, chr, pvalue)
    def rowInfo(self, v, tissueNr):
        inds = templ["dataInds"]
        return ([(int(v[inds["baitStart"]]), int(v[inds["baitEnd"]])), (int(v[inds["oeStart"]]), int(v[inds["oeEnd"]])), self.chInd["chr" + str(v[inds["baitChr"]])], float(v[inds["firstPvalue"] + tissueNr])])

    # Creates sorted lists of segments and links
    def initData(self):
        # Read input file
        df = pd.read_csv(self.inputPath, sep='\t')
        fdata = df.values
        for v in fdata:
            if str(v[templ["dataInds"]["baitChr"]]) == str(v[templ["dataInds"]["oeChr"]]) and "chr" + str(v[templ["dataInds"]["baitChr"]]) in self.chInd:
                # For each tissue with a high enough pvalue add both ends to segment list, and add tissue to linkTissues of current link
                for ind, tissue in enumerate(templ["tissueIDs"]):
                    if v[templ["dataInds"]["firstPvalue"] + ind] >= self.pvalue:
                        rv = self.rowInfo(v, ind)
                        self.chrs[rv[2]].segments.add(rv[0])
                        self.chrs[rv[2]].segments.add(rv[1])
                        self.chrs[rv[2]].setLinkInfo(rv, tissue)
        print("### processed:", self.inputPath)
        del df
        del fdata
        # Sort segments and assign IDs to them
        for c in self.chrNames:
            chr = self.chrs[self.chInd[c]]
            chr.segments = sorted(list(chr.segments))
            for s in chr.segments:
                if s in chr.segmInd:
                    print("!!!Segment", s, "already has id!!!")
                chr.indSegm(s)

    # Creates a dictionary of main data that should be written to json
    def finalData(self):
        rez = {"pvalue": self.pvalue, "chrNames": self.chrNames,
               "tissueIDs": self.tissueIDs, "tissueBits": self.tissueBits}
        chrs = {}
        for c in self.chrNames:
            chr = self.chrs[self.chInd[c]]
            # Create a list of links
            links = []
            for i in chr.linkTissues.keys():
                links.append(chr.getLinkRow(i))
            links.sort()
            segments = chr.segments
            # Check if order of segments matches data in segmInd dictionary
            for i in range(len(segments)):
                if i != chr.segmInd[segments[i]]:
                    print("### wrong order of segments", c, i, segments[i])
            chrs[c] = {"links": links, "segments": segments}
        rez["chrValues"] = chrs
        return rez


data = DataHiC(templ)
fdata = data.finalData()
if not os.path.exists(templ["outputPath"]):
    os.makedirs(templ["outputPath"])
outputPath = templ["outputPath"] + \
    "/data-bloodCellHiC-pv-" + str(data.pvalue) + "-fin.json"
with open(outputPath, 'w') as outfile:
    json.dump(fdata, outfile)
    print("DONE save ", outputPath)
