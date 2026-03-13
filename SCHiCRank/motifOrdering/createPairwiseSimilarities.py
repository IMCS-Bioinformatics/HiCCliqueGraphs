import pickle
import os
from itertools import combinations
import csv



def callPairwiseSimilarites(filename, work_dir=""):
    with open(filename, "rb") as f:
        data = pickle.load(f) #Read clique data
    iii=0
    typ = data["type"]

    motifNames = data["clique_cells"].keys()
    motifLengths = ["alllengths", "long"]
    cellIDs = data["cell_IDs"]

    for motifName in motifNames:
        for motifLength in motifLengths:

            resultDir = os.path.join(work_dir, f"pairwiseSimilarities/{motifName}-{motifLength}/")
            os.makedirs(resultDir, exist_ok=True)
            resultFn = os.path.join(resultDir, f'pairwiseSimilarities-{typ}-{data["chr"]}-{data["resolution"]}-{motifName}-{motifLength}.csv')
            print(f"Processing {resultFn}")
            cellPairFrequencies = dict()
            for clique, listOfCells in data["clique_cells"][motifName].items():
                if len(listOfCells) <=1:
                    continue

                #Check clique length
                if motifLength == "long":
                    if clique[-1]-clique[0]<2000000:
                        continue #Do not process this clique because it is considered short

                for cell1, cell2 in combinations(listOfCells, 2):
                    if cell1 > cell2:
                        cell1, cell2 = cell2, cell1
                    pair = (cell1, cell2)
                    if pair not in cellPairFrequencies:
                        cellPairFrequencies[pair] = 0
                    cellPairFrequencies[pair] += 1

            # Sort the pairs by frequency in descending order
            sorted_pairs = sorted(cellPairFrequencies.items(), key=lambda x: x[1], reverse=True)

            # Save the pairwise frequencies as csv
            with open(resultFn, "w", newline='') as csvfile:
                fieldnames = ["Item 1", "Item 2", "Frequency", "cell1_name", "cell2_name", "cell1_phase", "cell2_phase", "same?", "type"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for pair, frequency in sorted_pairs:
                    writer.writerow({
                        "Item 1": pair[0],
                        "Item 2": pair[1],
                        "Frequency": frequency,
                        "cell1_name": data["index_to_name"][pair[0]],
                        "cell2_name": data["index_to_name"][pair[1]],
                        "cell1_phase": data["index_to_type"][pair[0]],
                        "cell2_phase": data["index_to_type"][pair[1]],
                        "same?": data["index_to_type"][pair[0]] == data["index_to_type"][pair[1]],
                        "type": f'{data["index_to_type"][pair[0]]}+{data["index_to_type"][pair[1]]}'
                    })




if __name__ == "__main__":

    callPairwiseSimilarites("base100k-chr18-100000-cliques.pkl")
