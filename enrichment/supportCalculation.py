from universal import UniversalDS, ChrData
from topologicalFeatures import *

templ = {
    "outputPath": "./transformedData/CliqueSegmentStateCountNew/",
    "minIntersection": 0.2,
}

blood_data = {
    "dataSet": "BloodCellPCHi-C",
    "dataPath": "./transformedData/BloodCellPCHi-C/Universal/",
    "pvalues": [5],
    "outputPath": "./transformedData/BloodCellPCHi-C/CliqueBaseSegments/",
}

normalHiC_data = {
    "dataSet": "NormalHi-C",
    "dataPath": "./transformedData/NormalHi-C/Universal/",
    "pvalues": [10],
    "outputPath": "./transformedData/NormalHi-C/CliqueBaseSegments/",
}

pcHiC_data = {
    "dataSet": "PCHi-C",
    "dataPath": "./transformedData/PCHi-C/Universal/",
    "pvalues": [0.7],
    "outputPath": "./transformedData/PCHi-C/CliqueBaseSegments/",
}

data_sets = [blood_data, pcHiC_data, normalHiC_data]

for ds in data_sets:
    for pv in ds["pvalues"]:
        fn = ds["dataPath"] + "/data-pvalue-" + str(pv) + "-fin.json"
        outputPath = ds["outputPath"]
        U = UniversalDS(fn)

        baseNsegments = {}
        bases = [5, 15, 30, 50]
        for base in bases:
            baseNsegments[base] = {}

        for ch in U.chrs:
            chData = ChrData(U, ch=ch, minLinkTissueCount=1)
            Bor = Bases(chData)
            Bor.setBasesOr()
            for base in bases:
                baseNsegments[base][ch] = {}
                baseN = Bor.getDegreeNBases(base)
                for tis in baseN:
                    segments = set()
                    for b in baseN[tis]:
                        segments.add(b[0])
                        segments.add(b[1])
                    baseNsegments[base][ch][tis] = sorted(list(segments))
        with open(outputPath + "baseN-segments-OR-pvalue-" + str(pv) + ".json", "w") as outfile:
            json.dump(baseNsegments, outfile)
