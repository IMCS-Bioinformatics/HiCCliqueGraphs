import pandas as pd
import os
import json

templ = {
    "dataPath": "./sourceData/pcHi-C(hg19)/",
    "outputPath": "./transformedData/PCHi-C",
    "tissueIDs": ['AD2', 'AO', 'BL1', 'CM', 'EG2', 'FT2', 'GA', 'GM', 'H1', 'HCmerge', 'IMR90', 'LG', 'LI11', 'LV', 'ME', 'MSC', 'NPC', 'OV2', 'PA', 'PO3', 'RA3', 'RV', 'SB', 'SG1', 'SX', 'TB', 'TH1', 'X5628FC'],
    "dataInds": {"pvalue": 11, "from": 1, "to": 2},
    "pvalue": 0.7,
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
        s1 = (int(v[0].split("-")[0]), int(v[0].split("-")[1]))
        s2 = (int(v[1].split("-")[0]), int(v[1].split("-")[1]))
        # Sort link segments if necessary
        if s1 > s2:
            s1, s2 = s2, s1
        link = (s1, s2)
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
        self.tissueIDs = templ["tissueIDs"]
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
        for i in range(len(self.tissueIDs)):
            self.tissueInd[bitNames[i]] = i
        self.chrs = [DataChr(self) for i in range(len(self.chrNames))]
        self.initData()

    # Returns needed data from a row in PCHi-C data (one end segment, other end segment, chr, pvalue)
    def rowInfo(self, v):
        inds = templ["dataInds"]
        return ([v[inds["from"]].split(":")[1], (v[inds["to"]].split(":"))[1], self.chInd[v[inds["from"]].split(":")[0]], v[inds["pvalue"]]])

    # Creates sorted lists of segments and links
    def initData(self):
        # Read all zip files in input directory
        arr = os.listdir(self.inputPath)
        zf = [z for z in arr if '.txt.zip' in z]
        for fn in zf:
            df = pd.read_csv(self.inputPath + "/" + fn,
                             sep='\t', compression='zip')
            fdata = df.values
            tissue = (fn.split("."))[0]
            # For each row (interaction) with a high enough pvalue add both ends to segment list, and add tissue to linkTissues of given link
            for v in fdata:
                if v[templ["dataInds"]["pvalue"]] >= self.pvalue:
                    rv = self.rowInfo(v)
                    self.chrs[rv[2]].segments.add(
                        (int(rv[0].split("-")[0]), int(rv[0].split("-")[1])))
                    self.chrs[rv[2]].segments.add(
                        (int(rv[1].split("-")[0]), int(rv[1].split("-")[1])))
                    self.chrs[rv[2]].setLinkInfo(rv, tissue)
            print("### processed:", fn, tissue)
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
    "/data-pcHiC-pv-" + str(data.pvalue) + "-fin.json"
with open(outputPath, 'w') as outfile:
    json.dump(fdata, outfile)
    print("DONE save ", outputPath)
