import pickle
import networkx as nx
import os


def createCliquePickles(
                        baseFn: str,
                        resFn: str,
                        ):
    
    with open(baseFn, 'rb') as f:
        data = pickle.load(f) #Links without calculated cliques


    resObj = {
        "chr": data["chr"],
        "resolution": data["resolution"],
        "type": data["type"]+"_cliques",
        "index_to_name": data["index_to_name"],
        "index_to_type": data["index_to_type"],
        "cell_IDs": data["cell_IDs"],
        "cell_cliques": {f"K{N}":{cellID: set() for cellID in data["cell_IDs"]} for N in [3,4,5,6,7,8]}, #For each cell, a set of cliques this cell has
        "clique_cells": {f"K{N}":{} for N in [3,4,5,6,7,8]}, #For each clique, set of cells where this clique is found
    }

    skipped = 0
    for cellID in data["cell_IDs"]:
        if cellID not in data["cell_links"]:
            skipped += 1
            continue
        if (cellID + 1) % 100 == 0:
            print(f"  Cell {cellID + 1}/{len(data['cell_IDs'])}")
        G = nx.Graph()
        G.add_edges_from(data["cell_links"][cellID])

        cliques_up_to_8 = [clique for clique in nx.find_cliques(G) if len(clique) <= 8]
        for cliqueSize in [3,4,5,6,7,8]:
            KN = f"K{cliqueSize}"
            cliques = [tuple(sorted(clique)) for clique in cliques_up_to_8 if len(clique) == cliqueSize]
            resObj["cell_cliques"][KN][cellID] = set(cliques)
            for clique in cliques:
                if clique not in resObj["clique_cells"][KN]:
                    resObj["clique_cells"][KN][clique] = set()
                resObj["clique_cells"][KN][clique].add(cellID)
    if skipped > 0:
        print(f"  Warning: Skipped {skipped} cells with no links in {data['chr']}")

    with open(resFn, 'wb') as f:
        pickle.dump(resObj, f)
    print("Saved cliques to", resFn)
    return resObj



if __name__ == "__main__":
    # Example usage

    for ch in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"]:
        baseFn = f"base60k-{ch}-60000.pkl"
        resFn = f"base60k-{ch}-60000-cliques.pkl"
        createCliquePickles(baseFn, resFn)
    
