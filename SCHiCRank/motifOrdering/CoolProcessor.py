import cooler
import networkx as nx
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

import random 

def extract_filename(file_path):
    # Get the file name from the full path and remove the .cool extension
    file_name = os.path.basename(file_path.strip("'"))
    return os.path.splitext(file_name)[0]  # Remove the extension

def aggregate_dicts_to_csv(dict_list, output_file):
    """
    Aggregates a list of dictionaries with the same keys into a CSV file.

    Parameters:
    dict_list (list): List of dictionaries to be aggregated.
    output_file (str): The name of the output CSV file.

    Example:
    dict_list = [
        {'A': 1, 'B': 2, 'C': 3},
        {'A': 4, 'B': 5, 'C': 6},
        {'A': 7, 'B': 8, 'C': 9}
    ]
    aggregate_dicts_to_csv(dict_list, "output.csv")
    """
    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(dict_list)
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    
    print(f"Data has been saved to {output_file}")

class CoolProcessor:
    """
    CoolProcessor is a class for handling and processing Hi-C data stored in .cool files using the cooler library.
    Attributes:
        c (cooler.Cooler): The Cooler object representing the .cool file.
        originalPath (str): Path to the original .cool file.
        cellName (str): Name of the cell or sample.
        chromNames (set): Set of chromosome names present in the .cool file.
        coolMatrix: Matrix representation of the Hi-C data (unbalanced, sparse).
        bins (DataFrame): DataFrame of bin information from the .cool file.
        chrInteractions (dict): Dictionary mapping chromosome names to sets of interactions (tuples of bin indices and counts).
        translateBinToLocus (dict): Dictionary mapping chromosome names to dictionaries translating bin indices to genomic loci.
        allChrJsonInteractions (dict): Dictionary mapping chromosome names to lists of interactions with loci and counts.
    Methods:
        __init__(self, pathToCoolFile, cellName=None):
            Initializes the CoolProcessor with a path to a .cool file or a Cooler object.
        __setInteractions(self, ch):
            Internal method. Loads and caches all interactions for a given chromosome, and sets up bin-to-locus translation.
        getInteractions(self, ch):
            Returns a set of interactions (tuples of bin indices and counts) for the specified chromosome.
        hasDataOnChr(self, ch):
            Checks if the specified chromosome is present in the .cool file.
        reduceResolution(self, k, fn="tmp.cool"):
            Reduces the resolution of the .cool file by a coarsening factor k and returns a new CoolProcessor for the reduced file.
        __setAllInteractionsWithLoci(self, ch):
            Internal method. Loads and caches all interactions for a chromosome, translating bin indices to genomic loci and summing counts.
        getAllInteractionsWithLoci(self, ch):
            Returns a list of interactions for the specified chromosome, with bin indices translated to genomic loci and counts summed.
    """
    def __init__(self, pathToCoolFile, cellName=None):
        if type(pathToCoolFile)==str:
            self.c = cooler.Cooler(pathToCoolFile)
            self.originalPath = pathToCoolFile
        else:
            self.c=pathToCoolFile
            self.originalPath = ""
        if cellName is None:
            cellName = extract_filename(pathToCoolFile)
        self.cellName = cellName
        
        self.chromNames = set(self.c.chromnames)
        self.coolMatrix = self.c.matrix(balance=False, sparse=True)
        self.bins = self.c.bins()[:]

        self.chrInteractions = {ch: None for ch in self.chromNames}
        self.translateBinToLocus = {ch: None for ch in self.chromNames}

        self.allChrJsonInteractions = {ch: None for ch in self.chromNames}

        print("CoolProcessor constructed")
    

    def __setInteractions(self, ch):
        pixels = self.coolMatrix.fetch(ch)
        rows, cols, counts = pixels.row, pixels.col, pixels.data

        # Directly create integer tuples without redundant conversions
        interactions = set(zip(rows.astype(int), cols.astype(int), counts.astype(int)))
        self.chrInteractions[ch] = interactions

        # Pre-fetch bin starts as numpy array for efficient indexing
        bin_starts = self.bins['start'].to_numpy().astype(int)

        # Set translation dict from index to locus
        translate_dict = {}
        for A, B, _ in interactions:
            translate_dict[A] = bin_starts[A]
            translate_dict[B] = bin_starts[B]

        self.translateBinToLocus[ch] = translate_dict

        return translate_dict

    def getInteractions(self, ch):
        if ch not in self.chromNames:
            return set()
        if self.chrInteractions[ch] is None:
            self.__setInteractions(ch)
        return self.chrInteractions[ch]
    
    def hasDataOnChr(self, ch):
        return ch in self.chromNames
    
    def reduceResolution(self, k, fn="tmp.cool"):
        #K- coarsen factor, int. e.g. 5 to make bin size 5 times bigger
        if k==1:
            return self
        
        coarsen_factor = k
        chunksize = 10000000  # You can adjust this based on your file size and memory

        # Create new cooler file with reduced resolution
        cooler.coarsen_cooler(self.originalPath, fn, factor=coarsen_factor, chunksize=chunksize)
        return CoolProcessor(fn, cellName=self.cellName)


    def __setAllInteractionsWithLoci(self, ch):
        #Links are sorted after this
        links = self.getInteractions(ch) #set of tuples (A,B,Count) where A and B are indeces
        lociLinks = [[self.translateBinToLocus[ch][link[0]], self.translateBinToLocus[ch][link[1]], link[2]] for link in links]
        lociLinks = [(min(A,B), max(A,B), count) for [A,B,count] in lociLinks]
        lociLinks = list(set(lociLinks))
        lociLinks = sorted(lociLinks, key=lambda x:(x[0], x[1]))

        linkToTuples = {}
        for (A,B,count) in lociLinks:
            if (A,B) not in linkToTuples:
                linkToTuples[(A,B)] = []
            linkToTuples[(A,B)].append(count)
        lociLinks = [[A,B,sum(counts)] for (A,B), counts in linkToTuples.items()]
        lociLinks = sorted(lociLinks, key=lambda x:(x[0],x[1]))

        self.allChrJsonInteractions[ch] = lociLinks
        return lociLinks
    
    def getAllInteractionsWithLoci(self, ch):
        #Returns list of links [aaa, bbb, count] where aaa and bbb are loci
        if not self.allChrJsonInteractions[ch]:
            self.__setAllInteractionsWithLoci(ch)
        return self.allChrJsonInteractions[ch]
    