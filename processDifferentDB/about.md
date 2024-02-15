This directory contains scripts to collect C3 segment data for Hi-C graphs from Encode and ensembl databases.

This directory with all data files (that take much storage) is available at http://susurs.mii.lu.lv/HiCData/processDifferentDB.zip

Different Hi-C graphs for different tissue types are taken, cliques of size 3 (C3) are found in each graph, and segments that are a part of any C3 are analyzed, by finding data about their genomic loci in different databeses, such as encode, ensembl.

Script that generates overview files for 2 Hi-C datasets that were used in our studies - [Tissue PCHi-C by Jung et al](https://doi.org/10.1038/s41588-019-0494-8) and [Tissue Hi-C from 3DIV](https://doi.org/10.1093/nar/gkaa1078) - is found in [processEncodeEnsembl.py](./processEncodeEnsembl.py).

A sample result is also included [here](./graphFilesTissuePCHiC/rez/overview-Fantoms-GTEx-Segments-tissuePCHiC-chrX-v0205.csv). It shows information about every segment.
- isInClique: true if this node is a part of at least one C3 in the current tissue type graph.
Columns about FANTOM5 array of genome-wide transcription profiles (Forrest et al. (2014))
- fantomDistance: distance in base pairs to the closest FANTOM5 segment.
- fantomExpression: the expression value of the closest FANTOM5.
- fantomCount: numebr of FANTOM5 that were the closest to the given segment. If the value is greater than one and fantomDistance=0, then the current segment overlaps with more than one FANTOM5 segment, and the fantomExpression value is the max expression of all fantom5 segments.

GTExExpression, GTExDistance, GTExCount - same as for FANTOM5, but using Roadmap Epigenomics per Kundaje et al. (2015).
- RNAPolCount: 0 if there is no RNA Polymerase in the current segment. RNA polymerase II sites in our cliques are obtained using TF ChIP-seq data from ENCODE (see Luo et al. (2020))

The Hi-C data as well as the files from different databases for these datasets are too large to be included in Github, therefore this whole folder with all data files is available [here](http://susurs.mii.lu.lv/HiCData/processDifferentDB.zip).

The data for graph segments is taken from:
- FANTOM5 array of genome-wide transcription profiles (Forrest et al. (2014))
- Roadmap Epigenomics per Kundaje et al. (2015)
- TF ChIP-seq data from ENCODE (see Luo et al. (2020))


