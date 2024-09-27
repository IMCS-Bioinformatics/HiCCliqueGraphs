# ENCODE cCRE clique expression analysis
This directory contains scripts to analyze gene GTEx and FANTOM5 profile expressions in various types of cliques with cCRE annotations.

## Example datasets used
We calculated gene and FANTOM5 expressions for different clique types for 2 datasets: [Tissue pcHi-C by Jung et al](https://doi.org/10.1038/s41588-019-0494-8) and [Tissue Hi-C from 3DIV](https://doi.org/10.1093/nar/gkaa1078)

FANTOM5 array of genome-wide transcription profiles [Forrest et al. (2014)](https://www.nature.com/articles/nature13182)

Gene tissue expressions from GTEx [Lonsdale et al. (2013)](https://doi.org/10.1038/ng.2653)

Chromatin cCRE annotations from [ENCODE](https://doi.org/10.1093/nar/gkz1062)

## Overview of the scripts in this directory
There are 3 scripts here. They should be executed in this order, because the results of each script are used as input data for the next script

1) [prepareForcCRE.py](./prepareForcCRE.py) assigns cCRE annotations, gene and FANTOM5 expression profiles to Hi-C network nodes;
2) [finalProcess_cCRE.py](./finalProcess_cCRE.py) creates a csv file where gene and FANTOM5 profile expression statistics (mean, median, etc.) are calculated for different clique types;
3) [create_cCRE_plots.R](./create_cCRE_plots.R) uses the generated csv file and generates images that show differences in expressions for different clique types.

### Usage

Please clone the repository locally and run python modules from the [top directory of the repository](../)

Both python modules in the current directory contain templates at the start of the code that set all necessary variables for each of the example datasets. User should select the dataset by uncommenting the line for the desired dataset.

#### Assigning properties to graph nodes

The [prepareForcCRE.py](./prepareForcCRE.py) module requires Hi-C dataset in a preprocessed format and data files with ENCODE cCRE annotations, gene and FANTOM5 expression profile locations and their expressions. These datasets are too large to be pushed to GitHub, and they are available in our [dedicated page](http://susurs.mii.lu.lv/HiCData/). For demo purposes, the output of this module for 2 chormosomes is available in [results](./results/data-pvalue-0.7-fin-digraph-GRCh38-mini-wEncodeBits-wExpressions-wC3.json) and it can be used by other modules.

To run [prepareForcCRE.py](./prepareForcCRE.py) module that assigns cCRE annotations, gene and FANTOM5 profiles to Hi-C network nodes, please download:
1) Hi-C datasets in a preprocessed format for [tissue pcHi-C dataset](http://susurs.mii.lu.lv/HiCData/data-pvalue-0.7-fin-digraph-GRCh38.json) and [tissue pcHi-C dataset](http://susurs.mii.lu.lv/HiCData/data-pvalue-10-fin-GRCh38.json) and place them in [data directory](./data/);
2) Datafiles with cCRE annotations, gene and FANTOM5 expression profiles from [here](http://susurs.mii.lu.lv/HiCData/processDifferentDB.zip), then extract the folder and replace the [processDifferentDB](../processDifferentDB/) with the extracted folder. The downloaded folder contains all necessary datafiles.

The results will be saved in [results directory](./results/). As an example, [result](./results/) for processing chr5 and chr10 of tissue pcHi-C dataset are shown there and can be used by the other module.

#### Calculating expression statistics of different clique types
The [finalProcess_cCRE.py](./finalProcess_cCRE.py) module calculates statistics of gene and FANTOM5 profile expressions in different node subsets in different graphs.

To run, please uncomment the desired dataset template. A third (default) template is available that uses a [sample result](./results/data-pvalue-0.7-fin-digraph-GRCh38-mini-wEncodeBits-wExpressions-wC3.json) with 2 chromosomes and has a smaller file size. To obtain results for all chromosomes, please follow instructions on running [prepareForcCRE.py](./prepareForcCRE.py) with full datasets.

The result is a csv file for each dataset: [for tissue pcHi-C dataset](./results/stats-encodeCCRE-tisPCHiC.csv) and [for tissue Hi-C dataset](./results/stats-encodeCCRE-3DIV.csv)

#### Plotting images
To plot some aspects of the results, i.e. to see relationships between expressions of different clique types, use [R script](./create_cCRE_plots.R)

## Description of cliques that are analyzed

First, cCRE annotations are added to Hi-C network nodes. Each node represents a chromatin segment. There are 4 possible cCRE annotations: promoter (P), distal enhancer (E), CTCF (C) and DNase-only (D).
Next, genes and FANTOM5 profiles are added to graph nodes. The same principle of assigning biological properties to nodes is used - if node segment and property segment overlap by at least 20% of the shortest interval, a property is assigned to a node. An exception is FANTOM5 profile assignment, where it needs to be located not further than 2000bp away from a node.

We use GTEx and FANOTM5 data to obtain information about different gene and FANTOM5 profile expressions in different tissue types. 

In the Hi-C network for one tissue type and for one chromosome (we ignore interchromosomal interactions), we find different topological elements, e.g. cliques of size 3, 4 or 5, or cliques with one edge missing. For the topological elements we found, we analyze gene or FANTOM5 expressions: we collect a set of all genes that were assigned to the nodes in the topological elements, and we calculate statistics (mean, median, quartiles) of the expressions.

Alternatively, we can set stricter requirements for topological elements based on cCRE annotations in the nodes. For example, we may consider a set of topological elements (e.g. cliques of size 4) where all nodes must have distal enhancer annotation, or all nodes must have specific cCRE annotation combinations. We then similarly calculate expression statistics for the gene and FANTOM5 expressions for a specific subset of topological elements, and then compare the expressions in different topological elements.



