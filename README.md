# HiCCliqueGraphs
This is a repository for **Clique-based topological characterization of Hi-C interaction graphs.**

Scripts for calculating and analyzing different topological elements of Hi-C graphs are included here.

For randomization that produces topologically similar graphs, see [randomization](./randomization/fastRandomizer)

## Usage
Clone this repository to your local computer to use it.
Some scripts include hard-coded relative paths to sample data files. Please run .py scripts from the root directory of the repository.

### Dependencies
- python 3.9
- networkx 2.8.4
- goatools 1.2.4
- intervaltree 3.1.0
- pyliftover 0.4

These versions were used by us. Other versions are not guaranteed to work.

## Data availability
Sample datasets with data on two chromosomes are included in [sampleData](./sampleData) and all scripts in this repository use them as an example.
Full datasets that can be used with our scripts can be found [here](http://susurs.mii.lu.lv/HiCData/)

## Dataset comparison and C3 and S(k) calculation
Scripts to compare datasets and calculate cliques with 3 vertices (C3) and supports S(k) are located in [topological elements](./topologicalElements)


## Gene ontology enrichment analysis 
Scripts to perform GOEA are found in [GOEA](./GOEA) directory.

## Enrichment in chromatin annotations
Scripts to look for enrichment in certain chromatin annotations is found in [enrichment](./enrichment) directory.

## Randomization
In order to validate our results, we created graphs that are "similar" to original Hi-C interaction graphs. The randomized graphs have the following properties:
- The intersection of original edge set and the edge set of the randomized graph is less than some parameter;
- The number of edges is the same;
- Node degrees in both graphs are precizely the same;
- Distribution of chromatin interaction lengths in both graphs is as similar as possible. 

### Randomized graph generation
Script to generate a randomized graph that is "similar" to a given Hi-C graph can be found in [graph randomization](./randomization/fastRandomizer/) and is described in [our JBCB paper](10.1142/S0219720024400018)


## Images
Images used in [our ISBRA 2023 paper](https://doi.org/10.1007/978-981-99-7074-2_38) can be generated using scripts in [image generation](./images)

## Analysis of segments using Encode and Ensembl
Scripts to analyze the graph segments using additional biological data can be found in [processDifferentDB](./processDifferentDB)

## Analysis of gene and FANTOM5 expressions in cCRE cliques
Scripts for our ICCBB 2024 submission where we explore expressions of genes and FANTOM5 profiles in different topological elements can be found in [encodeCcreCliqueExpressions](./encodeCcreCliqueExpressions)

## Funding
This work is supported by the Latvian Council of Science project lzp-2021/1-0236
