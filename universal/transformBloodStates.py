import os
import csv

templ = {
    "dataPath": "./transformedData/",
    "statesPath": "./sourceData/chromatin states/RegBuild_BLUEPRINT_annotations/RegBuild_BLUEPRINT_annotations.txt",
    "outputPath": "./transformedData/RegBuild_BLUEPRINT_annotations_transformed/",
    "dataSet": "BloodCellPCHiC",
    "pvalue": 5,
    "minIntersection": 0.2
}

column_numbers = {
    "bait": {"chr": 1, "start": 2, "end": 3, "states": 11},
    "oe": {"chr": 5, "start": 6, "end": 7, "states": 13}
}
file = open(templ["statesPath"], "r")

stateFile = [[int(val) if val.isdigit() and i not in [1, 5] else val for i, val in enumerate(line.split('\t'))]
             for line in file]
# Only keep significant interactions in the same chromosome, and filter out chrY
stateFile = [line for line in stateFile if line[10]
             == 1 and line[1] == line[5] and line[1] != 'Y']

segmentStates = {}
chrs = set()
for line in stateFile:
    chr = line[1]
    chrs.add(chr)
    tissue = line[9]
    segmentStates[tissue] = segmentStates.get(tissue, {})
    segmentStates[tissue][chr] = segmentStates[tissue].get(chr, {})
    for interEnd in column_numbers:
        ind = column_numbers[interEnd]
        if line[ind["states"]] != '.':
            segm = tuple([line[ind["start"]], line[ind["end"]]])
            if segm not in segmentStates[tissue][chr]:
                segmentStates[tissue][chr][segm] = set()
                states = [state.split('_')[0]
                          for state in line[ind["states"]].split(',')]
                segmentStates[tissue][chr][segm].update(states)

chrs = sorted(list(chrs))
for tissue in segmentStates:
    transformedStates = []
    for chr in chrs:
        print(chr)
        chrTransformedStates = []
        for segm in segmentStates[tissue][chr]:
            for state in segmentStates[tissue][chr][segm]:
                chrTransformedStates.append(
                    ["chr" + chr, segm[0], segm[1], state])
        chrTransformedStates.sort()
        transformedStates.extend(chrTransformedStates)
    with open(f"{templ['outputPath']}/{tissue}_states.txt", 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(transformedStates)
