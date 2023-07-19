import itertools
import math
import os
import json
import CReader
import pandas as pd
import ntpath

templ = {
    "outputPath": "./results/EnrichmentInAnnotiations/",
    "minIntersection": 0.2,
}

blood_data = {
    "dataSet": "BloodCellPCHi-C",
    "dataPath": "./transformedData/BloodCellPCHi-C/Universal/",
    "supportPath": "./transformedData/BloodCellPCHi-C/CliqueBaseSegments/",
    "statesPath": "./transformedData/RegBuild_BLUEPRINT_annotations_transformed/",
    "pvalues": [5],
    "idToTissue": {"Ery": "Ery", "Mac0": "Mac0", "Mac1": "Mac1", "Mac2": "Mac2", "MK": "MK", "Mon": "Mon", "nCD4": "nCD4", "nCD8": "nCD8", "Neu": "Neu"},
    "stateCategories": {"ctcf": "ctcf", "distal": "distal", "dnase": "dnase", "proximal": "proximal", "tfbs": "tfbs", "tss": "tss"}
}

pcHiC_data = {
    "dataSet": "PCHi-C",
    "dataPath": "./transformedData/PCHi-C/Universal/",
    "supportPath": "./transformedData/PCHi-C/CliqueBaseSegments/",
    "statesPath": "./sourceData/chromatin states/18states.all.mnemonics.bedFiles",
    "pvalues": [0.7],
    "idToTissue": {"E080": "AD2", "E065": "AO", "E079": "EG2", "E063": "FT2", "E094": "GA", "E116": "GM", "E003": "H1", "E071": "HCmerge", "E126": "IMR90", "E096": "LG", "E066": "LI11", "E095": "LV",
                   "E004": "ME", "E006": "MSC", "E007": "NPC", "E097": "OV2", "E098": "PA", "E100": "PO3", "E104": "RA3", "E105": "RV", "E109": "SB", "E106": "SG1", "E113": "SX", "E005": "TB", "E112": "TH1", "E073": "X5628FC"},
    "stateCategories": {'1_TssA': 'TSS', '2_TssFlnk': 'TSS', '3_TssFlnkU': 'TSS', '4_TssFlnkD': 'TSS', '5_Tx': 'TX', '6_TxWk': 'TX', '7_EnhG1': 'ENH', '8_EnhG2': 'ENH', '9_EnhA1': 'ENH',
                        '10_EnhA2': 'ENH', '11_EnhWk': 'ENH', '12_ZNF/Rpts': 'ZNF/RPTS', '13_Het': 'HET', '14_TssBiv': 'TSS_BIV', '15_EnhBiv': 'ENH_BIV', '16_ReprPC': 'PC', '17_ReprPCWk': 'PC', '18_Quies': 'QUI'}
}

normalHiC_data = {
    "dataSet": "NormalHi-C",
    "dataPath": "./transformedData/NormalHi-C/Universal/",
    "supportPath": "./transformedData/NormalHi-C/CliqueBaseSegments/",
    "statesPath":  "./sourceData/chromatin states/18states.all.mnemonics.bedFiles",
    "pvalues": [10],
    "idToTissue": {"E080": "AD", "E065": "AO", "E125": "ASCE", "E003": "H1", "E126": "F0", "E096": "LG", "E066": "LI",
                   "E004": "ME", "E006": "MS", "E007": "TP", "E097": "OV", "E098": "PA", "E100": "PO", "E109": "SB", "E113": "SX", "E112": "TH"},
    "stateCategories": {'1_TssA': 'TSS', '2_TssFlnk': 'TSS', '3_TssFlnkU': 'TSS', '4_TssFlnkD': 'TSS', '5_Tx': 'TX', '6_TxWk': 'TX', '7_EnhG1': 'ENH', '8_EnhG2': 'ENH', '9_EnhA1': 'ENH',
                        '10_EnhA2': 'ENH', '11_EnhWk': 'ENH', '12_ZNF/Rpts': 'ZNF/RPTS', '13_Het': 'HET', '14_TssBiv': 'TSS_BIV', '15_EnhBiv': 'ENH_BIV', '16_ReprPC': 'PC', '17_ReprPCWk': 'PC', '18_Quies': 'QUI'}
}

data_sets = [blood_data, pcHiC_data, normalHiC_data]

##########################
# Assumptions: each set of segments (base and state) are disjoint, therefore one parallel traversal is enough
# Base segments are sorted
# State files are grouped by chr, and sorted by segments
##########################


# Returns a list of paths for all files in a directory
def dirFileList(dir):
    filenames = [os.path.join(dir, fn) for fn in os.listdir(dir)]
    return filenames


# Finds links (interactions) present in all given tissues, and chromatin segments with at least one link present in all tissues
def getTissueBase(allBaseSegments, links, tissueBit):
    tissueBaseSegmentIDs = set()
    baseLinks = []
    for link in links:
        if link[2] & tissueBit == tissueBit:
            baseLinks.append(link[0:2])
            tissueBaseSegmentIDs.update(link[0:2])
    baseSegments = [allBaseSegments[id] for id in tissueBaseSegmentIDs]
    baseSegments.sort()
    return baseSegments, baseLinks


# Creates a dictionary {base segment: list of states}
# Assumes all segments are sorted by start
def getSegmentStates(baseSegments, stateSegments):
    assert baseSegments == sorted(baseSegments), "base segments not sorted"
    assert stateSegments == sorted(stateSegments), "state segments not sorted"
    segmentStates = {}
    stateSegmIter = 0
    currStateSegments = set()
    for baseSegm in baseSegments:
        segmentStates[tuple(baseSegm)] = set()
        while stateSegmIter < len(stateSegments) and stateSegments[stateSegmIter][1] <= baseSegm[1]:
            currStateSegments.add(tuple(stateSegments[stateSegmIter]))
            stateSegmIter += 1
        currStateSegmList = list(currStateSegments)
        for stateSegm in currStateSegmList:
            if stateSegm[2] < baseSegm[0]:
                currStateSegments.remove(stateSegm)
            else:
                if intersectsEnough(stateSegm[1: 3], baseSegm, templ["minIntersection"]):
                    state = data_set["stateCategories"][stateSegm[3]]
                    segmentStates[tuple(baseSegm)].add(state)

    return segmentStates


# Checks if intersection of two segments is at least minIntersection of smallest segment
def intersectsEnough(segm1, segm2, minIntersection):
    if not ((segm1[1] <= segm2[0]) or (segm2[1] <= segm1[0])):
        allPoints = [segm1[0], segm1[1], segm2[0], segm2[1]]
        allPoints.sort()
        if (allPoints[2] - allPoints[1]) / min(segm1[1] - segm1[0], segm2[1] - segm2[0]) >= minIntersection:
            return True
    return False


# Creates an adjacency list for the graph
def getAdj(links):
    adj = {}
    for link in links:
        for i in [0, 1]:
            if link[i] not in adj:
                adj[link[i]] = set()
            adj[link[i]].add(link[1-i])
    return adj


# Creates a list of segments that are in cliques (C3)
def getClique3SegmentIDs(segments, links):
    cliqueSegments = set()
    adj = getAdj(links)
    for link in links:
        if not (link[0] in cliqueSegments and link[1] in cliqueSegments):
            if len(adj[link[0]] & adj[link[1]]) > 0:
                cliqueSegments.add(link[0])
                cliqueSegments.add(link[1])
    return sorted(list(cliqueSegments))


# For each segment creates a set of states that are present in all tissues
def getSegmentStatesAnd(segmentStates, tissues):
    segmStatesAnd = {}
    for segm in segmentStates[tissues[0]]:
        for state in segmentStates[tissues[0]][segm]:
            stateOk = True
            for tissue in tissues:
                if segm not in segmentStates[tissue] or state not in segmentStates[tissue][segm]:
                    stateOk = False
                    break
            if stateOk:
                segmStatesAnd[segm] = segmStatesAnd.get(
                    segm, set()) | set([state])
    return segmStatesAnd


# Counts in how many segments a state is present in cliques/supports vs general
def getStateCounts(segmentStates, cliqueSegments, tissues):
    stateCountForLinks = {}
    stateCountForCliques = {}
    for segm in segmentStates:
        for state in segmentStates[segm]:
            stateCountForLinks[state] = stateCountForLinks.get(
                state, 0) + 1
            if segm in cliqueSegments:
                stateCountForCliques[state] = stateCountForCliques.get(
                    state, 0) + 1
    return stateCountForLinks, stateCountForCliques


# Adds a row to result data, including ratio of state count in cliques vs general
def updateResults(support, segmentCountForLinks, cliqueSegments, tissues):
    segmentCountForCliques = len(cliqueSegments)
    if support not in dataForCliques:
        dataForCliques[support] = []
        dataForCliquesNorm[support] = []
        dataForRatios[support] = []
    stateCounts = [stateCountForCliques.get(state, 0) for state in categories]
    dataForCliques[support].append(tissues +
                                   [chr, segmentCountForCliques] + stateCounts)

    stateCountsNorm = [stateCountForCliques.get(
        state, 0) / segmentCountForCliques if segmentCountForCliques > 0 else -0.01 for state in categories]
    dataForCliquesNorm[support].append(tissues +
                                       [chr, segmentCountForCliques] + stateCountsNorm)

    stateRatios = [(stateCountForCliques.get(state, 0) / segmentCountForCliques)
                   / (stateCountForLinks.get(state, 0) / segmentCountForLinks)
                   if stateCountForLinks.get(state, 0) > 0 and segmentCountForCliques > 0
                   else -0.01 for state in categories]
    dataForRatios[support].append(tissues +
                                  [chr, segmentCountForCliques / segmentCountForLinks if segmentCountForLinks > 0 else -0.01] + stateRatios)


# Gets tissue from filename (assuming it is right before the first underscore)
def getTissueIdFromFilename(dataSet, filename):
    id = ntpath.basename(filename).split("_")[0]
    tissue = data_set["idToTissue"].get(id, "")
    return id, tissue


# Writes results to spreadsheets with formatting
def writeResults(fn):
    supports = sorted(dataForCliques.keys(), key=lambda x: int(x))
    with pd.ExcelWriter(fn) as writer:
        percent_format = writer.book.add_format(
            {'num_format': '0%'})
        df = pd.DataFrame(dataForLinks, columns=columns)
        df.to_excel(writer, sheet_name='Links')
        df = pd.DataFrame(dataForLinksNorm, columns=columns)
        df.to_excel(writer, sheet_name='LinksNorm')
        writer.sheets['LinksNorm'].conditional_format(1, 3 + len(tissues), 1 + math.comb(len(allTissues), tissueCount) * 23, 2 + len(categories) + len(tissues),
                                                      {'type': '3_color_scale', 'minimum': 0, 'maximum': 1})
        writer.sheets['LinksNorm'].set_column(
            3 + len(tissues), 2 + len(categories) + len(tissues), None, percent_format)
        for support in supports:
            columns[tissueCount + 1] = "Total support segment count"
            df = pd.DataFrame(dataForCliques[support], columns=columns)
            df.to_excel(writer, sheet_name='Support' + support)
            df = pd.DataFrame(
                dataForCliquesNorm[support], columns=columns)
            columns[tissueCount + 1] = "Support to all segment ratio"
            df.to_excel(writer, sheet_name='Support' + support + 'norm')
            df3 = pd.DataFrame(dataForRatios[support], columns=columns)
            df3.to_excel(writer, sheet_name='Ratio' + support)
            writer.sheets['Support' + support + 'norm'].conditional_format(1, 3 + len(tissues), 1 + math.comb(len(allTissues), tissueCount) * 23, 2 + len(categories) + len(tissues),
                                                                           {'type': '3_color_scale', 'minimum': 0, 'maximum': 1})
            writer.sheets['Support' + support + 'norm'].set_column(
                3 + len(tissues), 2 + len(categories) + len(tissues), None, percent_format)

            writer.sheets['Ratio' + support].conditional_format(1, 3 + len(tissues), 1 + math.comb(len(allTissues), tissueCount) * 23, 2 + len(categories) + len(tissues),
                                                                {'type': '3_color_scale', 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num', 'min_value': 0, 'max_value': 2.5, 'mid_value': 1})
            writer.sheets['Ratio' + support].set_column(
                2 + len(tissues), 2 + len(categories) + len(tissues), None, percent_format)


for data_set in data_sets:
    categories = sorted(list(set(data_set["stateCategories"].values())))
    for pvalue in data_set["pvalues"]:
        baseData = CReader.readJson(
            data_set["dataPath"] + "data-pvalue-" + str(pvalue) + "-fin.json")
        segmentStates = {chr: {} for chr in baseData["chrNames"]}
        # Calculate for supports in tissue pairs and cliques (C3) in single tissues
        # for orAnd in ["OR", "AND", "C3"]:
        for orAnd in ["OR", "C3"]:
            # Support segments should already be calculated, but C3 segments will be calculated later
            if orAnd == "C3":
                supportNSegments = {"1": {chr: {}
                                          for chr in baseData["chrNames"]}}
            else:
                supportNpath = data_set["supportPath"] + orAnd + \
                    "_rezBaseIndexes2-pv" + str(pvalue) + ".json"
                supportNSegments = CReader.readJson(supportNpath)

            # Find states for segments
            filenames = dirFileList(data_set["statesPath"])
            allTissues = []
            for filename in filenames:
                id, tissue = getTissueIdFromFilename(data_set, filename)
                if id not in data_set["idToTissue"]:
                    continue
                allTissues.append(tissue)
                print("start", tissue)
                stateFile = [[int(val) if val.isdigit() else val for val in line.split()]
                             for line in open(filename, "r")]
                tissueBit = baseData["tissueBits"][tissue]
                for chr in baseData["chrValues"]:
                    stateSegments = [
                        line for line in stateFile if line[0] == chr]
                    links = baseData["chrValues"][chr]["links"]
                    segments = baseData["chrValues"][chr]["segments"]
                    tissueBaseSegments, links = getTissueBase(
                        segments, links, tissueBit)
                    segmentStates[chr][tissue] = getSegmentStates(
                        tissueBaseSegments, stateSegments)
                    if orAnd == "C3":
                        supportNSegments["1"][chr][tissue] = getClique3SegmentIDs(
                            segments, links)

            tissueCount = 1 if orAnd == "C3" else 2
            columns = ["Tissue"] * tissueCount + \
                ["Chr", "Total segment count"] + categories

            # Count states present in each base and support segment in all tissue combinations
            dataForLinks = []
            dataForLinksNorm = []
            dataForCliques = {}
            dataForCliquesNorm = {}
            dataForRatios = {}
            for chr in baseData["chrValues"]:
                segments = baseData["chrValues"][chr]["segments"]
                for tissues in itertools.combinations(allTissues, tissueCount):
                    segmStatesAnd = getSegmentStatesAnd(
                        segmentStates[chr], tissues)
                    for support in supportNSegments:
                        tissueBits = sum([baseData["tissueBits"][t]
                                         for t in tissues])
                        tissueSegments, links = getTissueBase(
                            segments, baseData["chrValues"][chr]["links"], tissueBits)
                        cliqueSegments = set([tuple(segm[0:2])
                                             for segm in tissueSegments])
                        for t in tissues:
                            cliqueSegments.intersection_update(set([tuple(segments[idx][0:2])
                                                                    for idx in supportNSegments[support][chr][t]]))
                        stateCountForLinks, stateCountForCliques = getStateCounts(
                            segmStatesAnd, cliqueSegments, tissues)
                        updateResults(
                            support, len(tissueSegments), cliqueSegments, list(tissues))
                    # Update link results
                    segmentCountForLinks = len(
                        tissueSegments)
                    stateCounts = [stateCountForLinks.get(
                        state, 0) for state in categories]
                    dataForLinks.append(
                        list(tissues) + [chr, segmentCountForLinks] + stateCounts)
                    # Add normalized data
                    stateCounts = [stateCountForLinks.get(
                        state, 0) / segmentCountForLinks if segmentCountForLinks > 0 else -0.01 for state in categories]
                    dataForLinksNorm.append(
                        list(tissues) + [chr, segmentCountForLinks] + stateCounts)

            # Write data to spreadsheet
            ds = data_set["dataSet"]
            if not os.path.exists(templ["outputPath"]):
                os.makedirs(templ["outputPath"])
            output_filename = f"{ds}-pv{pvalue}-stateSegmentCount-linksAndC3.xlsx" if orAnd == "C3" else f"{ds}-pv{pvalue}-stateSegmentCount-supportN-{orAnd}.xlsx"
            writeResults(templ["outputPath"] + output_filename)
