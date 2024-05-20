# ComponentVis pre-processing
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

## Flowchart
An image below shows the steps and data sets used during pre-processing to generate files that can be visualized using ComponentVis. The whole process is automated using bash scripts for the example data sets. 

![Flowchart for ComponentVis](https://github.com/IMCS-Bioinformatics/HiCCliqueGraphs/assets/119489036/80d2dd3b-dbb7-4372-9185-588cb0017a98)

## Using TADs calculated by other methods
It is possible to use TADs calculated by any method as long as they fit the following format. Boundaries of TADs should be stored in a JSON file and included in a ZIP file along with other pre-processed data. To include your own TAD files in the ZIP files, you can either place the TAD files in a directory "Processed_data/NormalHiC/TADs" (relative to bash script location) before running the pre-processing script or manually replace the included TAD files with your own after the ZIP files are created. The exact name of the TAD files does not matter; however, it must end with "Tad.json". Each file has to contain data for only one chromosome. The main JSON object should have two properties - "TADs" and "Real". "TADs" contain key-value pairs where each key is a tissue ID, and each value is an array containing all endpoints of TADs. The array must be sorted, and it must start with 0 and end with a number greater or equal to the length of the chromosome. "Real" contains key-value pairs where each key is a tissue ID, and each value is an array of zeros and ones. `data["TADs"][tissueID][i]` should be 1 if there is a TAD that starts at `data["Real"][tissueID][i]` and ends at `data["Real"][tissueID][i+1]`. Otherwise, it should be 0.

For example, if tissue with id "AD" had only two TADs `[40000, 480000]` and `[480000, 2645000]` and tissue with id "AO" also had only two TADs `[15000, 1360000]` and `[1715000, 1975000]`, then it could be represented as:
```json
{
  "TADs": {
    "AD": [0, 40000, 480000, 2645000, 500000000],
    "AO": [0, 15000, 1360000, 1715000, 1975000, 500000000]
  },
  "Real": {
    "AD": [0, 1, 1, 0],
    "AO": [0, 1, 0, 1, 0]
  }
}
```
