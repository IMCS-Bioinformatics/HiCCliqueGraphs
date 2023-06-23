# ISBRAprivate
This is a repository for our paper **Clique-based topological characterization of Hi-C interaction graphs**
Scripts for calculating and analyzing different topological elements of Hi-C graphs are included here.

## Usage
Clone this repository to your local computer to use it.

### Dependencies
- numpy
- matplotlib.pyplot
- networkx
- goatools
- intervaltree

## Data availability
Sample datasets with data on two chromosomes are included in [sampleData](./sampleData) and all scripts in this repository use them as an example.
Full datasets that can be used with our scripts can be found ....

## Dataset comparison and C3 and S(k) calculation
Scripts to compare datasets and calculate cliques with 3 vertices (C3) and supports S(k) are located in [topological elements](./topologicalElements)


## Gene ontology enrichment analysis 
Scripts to perform GOEA are found in [GOEA](./GOEA) directory.

## Enrichment in chromatin annotations

## Randomization
In order to validate our results, we created graphs that are "similar" to original Hi-C interaction graphs. The randomized graphs have the following properties:
- The intersection of original edge set and the edge set of the randomized graph is less than some parameter;
- The number of edges is the same;
- Node degrees in both graphs are precizely the same;
- Distribution of chromatin interaction lengths in both graphs is as similar as possible. 

### Randomized graph generation
Script to generate a randomized graph that is "similar" to a given Hi-C graph can be found in [graph randomization](./randomization/randomizer)
### Randomized graph analysis
Result validation by comparing real and generated graphs is found in [randomized graph analysis](./randomization/analysis)

## Images
Images used in our paper can be generated using scripts in [image generation](./images)

