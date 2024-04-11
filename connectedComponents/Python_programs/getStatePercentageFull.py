import os
import json
# import xlsxwriter
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--pvalue', type=float, default=6)
parser.add_argument('--dataPath', default=os.path.join(os.path.dirname(
    __file__), "../Processed_data/Normal_Hi-C/"))
args = parser.parse_args()

pvalue = args.pvalue
if not isinstance(pvalue, int) and pvalue.is_integer():
    pvalue = int(pvalue)

inputPath = args.dataPath + \
    "Interactions_with_states/Pvalue" + str(pvalue)
outputPath = args.dataPath + "State_percentage_stats/Pvalue" + str(pvalue)
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

countHas = {}
segmentsWithStates = {}
chrTissueSegmentCount = {}
allStates = set()
chrs = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
for chr in chrs:
    countHas[chr] = {}
    chrData = json.load(open(os.path.join(inputPath, chr + ".json")))
    tissues = chrData["tissueData"]["tissueBits"].keys()

    for tissue in tissues:
        countHas[chr][tissue] = {}
    tissueSegmentIDs = {}
    for tissue in tissues:
        tissueSegmentIDs[tissue] = set()
        segmentsWithStates[tissue] = {}
    edges = chrData["links"]
    nodes = chrData["segments"]
    tissueBits = chrData["tissueData"]["tissueBits"]
    for edge in edges:
        for tissue in tissueBits:
            if edge["bits"] & tissueBits[tissue] == tissueBits[tissue]:
                for segmentID in [edge["bin1"], edge["bin2"]]:
                    tissueSegmentIDs[tissue].add(segmentID)
                    node = nodes[segmentID]
                    if tissue in node["states"]:
                        for state in node["states"][tissue]:
                            if state not in segmentsWithStates[tissue]:
                                segmentsWithStates[tissue][state] = set()
                            segmentsWithStates[tissue][state].add(segmentID)
    chrTissueSegmentCount[chr] = {}
    for tissue in tissueBits:
        chrTissueSegmentCount[chr][tissue] = len(tissueSegmentIDs[tissue])
        for state in segmentsWithStates[tissue]:
            allStates.add(state)
            countHas[chr][tissue][state] = len(
                segmentsWithStates[tissue][state])


# Uncomment to also create a .xlsx file
# workbook = xlsxwriter.Workbook(os.path.join(
#     outputPath, 'statePercentageStats.xlsx'))
# worksheet = workbook.add_worksheet()
allStates = sorted(list(allStates))
# for i, state in enumerate(allStates):
#     worksheet.write(0, i+2, state)
fractions = {}
row = 1
for chr in chrs:
    fractions[chr] = {}
    for tissue in countHas[chr]:
        fractions[chr][tissue] = {}
#         worksheet.write(row, 0, chr)
#         worksheet.write(row, 1, tissue)
        for i, state in enumerate(allStates):
            if state not in countHas[chr][tissue]:
                countHas[chr][tissue][state] = 0
            fractions[chr][tissue][state] = countHas[chr][tissue][state] / \
                chrTissueSegmentCount[chr][tissue]
#             worksheet.write(
#                 row, i+2, (countHas[chr][tissue][state] / chrTissueSegmentCount[chr][tissue]))

#             worksheet.write(0, i+2, state)
        row += 1

# workbook.close()


with open(os.path.join(outputPath, 'statePercentageStats.json'), "w") as outputFile:
    json.dump(fractions, outputFile)
