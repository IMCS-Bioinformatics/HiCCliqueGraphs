# HiCCliqueGraphs
This directory contains all pre-processing scripts and programs for [ComponentVis](https://github.com/IMCS-Bioinformatics/ChromatinNetworkVisualisation/tree/main/Visualizations/ComponentVisualization).

[process-hic](./process-hic) includes code for transforming Hi-C data from its original format to a format that is used for other programs and adding ChromHMM annotations, genes and proteins.

Code for finding connected components is in [Find_components](./Find_components).

[Python_programs](./Python_programs) contains code for calculating additional data about the connected components and dataset in general, as well as combining multiple files for the same chromosome into ZIP files.

## Usage
Clone this directory to your local computer to use it.

downloadSourceData.py contains a script that creates the intended directories and downloads input data from our server.

Included bash scripts allow to more easily run all programs sequentially. The bash scripts should be run from this directory.

### Dependencies
The pre-processing has been tested using these versions of dependencies:
- Rust 1.67.1
- C++11 using g++ 9.2.0
- Python 3.12
