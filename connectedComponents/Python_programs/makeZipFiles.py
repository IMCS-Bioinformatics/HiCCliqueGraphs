from zipfile import ZipFile
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

inputJSON = args.dataPath + \
    "Interactions_with_states/Pvalue" + str(pvalue) + "/"
inputCSV = args.dataPath + "Components_extra/Pvalue" + str(pvalue) + "/"
inputTADs = args.dataPath + "TADs/"
inputStateStats = args.dataPath + \
    "State_percentage_stats/Pvalue" + str(pvalue) + "/"

outputZIP = args.dataPath + "ZIPs/" + "Pvalue" + str(pvalue) + "/"


chrs = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
chrs.reverse()
for chr in chrs:
    outputFilename = outputZIP + chr + ".zip"
    if not os.path.exists(outputZIP):
        os.makedirs(outputZIP)

    # Creating the ZIP file
    with ZipFile(outputFilename, "w") as outputFile:
        outputFile.write(inputJSON + chr + ".json", chr + ".json")
        outputFile.write(inputCSV + chr + ".csv", chr + ".csv")
        if os.path.exists(inputTADs + "RealRes-" + chr + "-SpectralTad.json"):
            outputFile.write(inputTADs + "RealRes-" + chr +
                             "-SpectralTad.json", "RealRes-" + chr + "-SpectralTad.json")
        if os.path.exists(inputStateStats + "statePercentageStats.json"):
            outputFile.write(
                inputStateStats + "statePercentageStats.json", "statePercentageStats.json")
