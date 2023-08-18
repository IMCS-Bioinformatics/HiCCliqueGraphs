# Fast randomization
This randomization takes a HiC graph as input and randomizes it while preserving most basic topological features - node degrees and link lengths.

Graph: nodes are chromatin segments and links are HiC contacts.

Link lengths are calculated in bin pairs between interacting segments.

## How to run?
The [main.py](./main.py) file contains a launcher for the randomization. It takes all parameters from a template file like [templateCool.json](./templateCool.json)

### template structure
At this moment, the randomization accepts 2 input formats: .cool file and link list.

#### Randomizing a link list
This is the preferred version.
If the Hi-C graph is given in a .csv file as a list of interactions like in [sampleOriginalLinks.csv](./sampleOriginalLinks.csv), where in each row one link is coded, delimited by commas, then the template should look like [templateLinks.json](./templateLinks.json).
- "useLinkList": true, #if false, a cool file will be expected
- "inputFn": "./sampleOriginalLinks.csv", #path to file with links
- "outputFn": "./myRandomizedLinks.csv",  #path to file where results will be stored
- "q": 0.7, #fraction of links to randomize in range (0,1)
- "runNaiveRandomization": false, #must be false most of the times. If true, a naive version of the randomziation will be used, where link lengths are ignored
- "coolRelated": {} #will be ignored if link list is given
 
#### Randomizing .cool file
In this case both input and output is a .cool file.
One chromsome is randomized, the others are copied without changes.
Usually there are many insignificant interactions in .cool files, whereas only significant interactions above some observed count treshold are chosen for the graph.
An option is provided (and suggested to use) to randomize only significant interactions, with count > some treshold. The resulting .cool file will still have all insignificant interactions (below treshold), but significant interactions will be randomized.

template.json structure (like in [templateCool.json](./templateCool.json)):
- "useLinkList": false, #a cool file is expected
- "inputFn": "./sampleCool.cool", #path to cool file
- "outputFn": "./sampleCool-randomized.cool", #path to result file
- "q": 0.7, #float in range (0,1) for fraction of links to swap
- "runNaiveRandomization": false,
- "coolRelated":{
    "chr": "15",
    "countTreshold": 4
} #chromosome and the count treshold. Note, in some .cool files chromosomes are named '1', '2',... while in others - 'chr1', 'chr2', ...


### Sample datasets in this directory
Some sample input files are provided for testing purposes.
- sampleCool.cool - [Schwarzer](Wibke Schwarzer et al. Two independent modes of chromatin organization revealed by cohesin removal.Nature, Nov 2017. 551(7678):51–56) dataset, with all chromosomes. Example usage in [templateCool](./templateCool.json)
- schwarzerChr15Sample.cool - chr15 only with filtered interactions having >4 observed count. Used in our abstract. Example usage in [templateCool2](./templateCool2.json)
- sampleOriginalLinks.csv - link list from schwarzerChr15Sample.cool. Example usage in [templateLinks](./templateLinks.json)





