import copy, json, math
import random
import numpy as np


def create_directory_if_not_exists(path):
    import os
    if not path: return
    if not os.path.exists(path):
        os.makedirs(path)

def createNoise(originalFileName, ch, NOISE_TO_ADD=0.2, RANDOM_LIGATION_NOISE=0.5, ITER=0, resDir="", DSname="tissue pcHiC"):
    """
    Creates noise within the data derived from an original file to simulate variations.

    Args:
        originalFileName (str): The name of the original file from which the data is taken.
        ch (str): The chromosome name, e.g. chr8
        NOISE_TO_ADD (float): The proportion of total noise to add to the data. 
        RANDOM_LIGATION_NOISE (float): The proportion of random ligation noise to add, in range [0,1]. The GENOMIC_EFFECT_NOISE will be 1-RANDOM_LIGATION_NOISE
        ITER (int): Iteration number, will be appended to the file name of the result
        resDir (str, optional): The directory path where the result files will be stored.
        DSname (str, optional): The name of the dataset, used for identification purposes.

    Returns:
        None. The function saves the results to the specified directory to a file with a name encapsulating all parameters that define the noise added to the dataset.

    """
    create_directory_if_not_exists(resDir)
    GENOMIC_EFFECT_NOISE = 1-RANDOM_LIGATION_NOISE
    noisyFn = f"{resDir}/noised-{DSname}-demo-{ch}-{NOISE_TO_ADD}-{RANDOM_LIGATION_NOISE}-{ITER}.json"
    with open(originalFileName, 'r') as f:
        U = json.load(f) #original file
    R = copy.deepcopy(U)

    oallLinks = U["chrValues"][ch]["links"] # Original all links
    oallSegments = U["chrValues"][ch]["segments"] # Original all segments

    
    midpoints = [(seg[1]+seg[0])//2 for seg in oallSegments]
    diffs = [midpoints[i+1]-midpoints[i] for i in range(len(midpoints)-1)]
    
    if min(diffs)>80000:
        minDiff=100000 
    else:
        minDiff = min(diffs)
        minDiff = 10**math.floor(math.log10(minDiff)) #either 10, 1000, 10000, ... 
        #This is to tackle the problem of varying bin sizes - all bins are approximated to one same size such that no 2 bins are located in the same location after approximation
    roundedSegments = [[seg[0],seg[1], int((seg[0]+seg[1])/(2*minDiff))*minDiff] for seg in oallSegments]
    iii=1
    thirds = [el[2] for el in roundedSegments]
    if len(set(thirds))!=len(roundedSegments):
        raise Exception("There are some segments that get rounded to the same value")

    segmentStartToIndex = {roundedSegments[i][2]:i for i in range(len(roundedSegments))}
    
    if "tissueBits" not in U: 
        tissueBits = U["tissueData"]["tissueBits"]
    else: 
        tissueBits = U["tissueBits"]
    randomizedLinkBits = {} #(A,B):bits where A and B are indeces of segments  
    for tis,tisBit in tissueBits.items():
        otissueLinks = [link for i,link in enumerate(oallLinks) if ((link[2]&tisBit)==tisBit)]
        olLengths = [abs(roundedSegments[B][2]-roundedSegments[A][2]) for [A,B,*pvs] in otissueLinks] 
        actualOLinks = [[roundedSegments[A][2],roundedSegments[B][2]] for [A,B,*pvs] in otissueLinks]
        noiseToAdd = int(NOISE_TO_ADD*len(otissueLinks)) #this many links to add
        rlnoiseToAdd = int(noiseToAdd*RANDOM_LIGATION_NOISE) #preserving node degrees
        gdnoiseToAdd = noiseToAdd-rlnoiseToAdd #preserving link lengths
        MAX_LINK_LENGTH = max([abs(el[1]-el[0]) for el in actualOLinks])

        #adding noise random ligation - node degrees
        RLNOISE = []
        nodeDegrees = {}
        for link in actualOLinks:
            if link[0] not in nodeDegrees:
                nodeDegrees[link[0]]=0
            nodeDegrees[link[0]]+=1
            if link[1] not in nodeDegrees:
                nodeDegrees[link[1]]=0
            nodeDegrees[link[1]]+=1
        items = list(nodeDegrees.keys()) #these are nodes
        weights = list(nodeDegrees.values()) #their relative weights - node degrees
        sampledItems = random.choices(items, weights=weights, k=rlnoiseToAdd)
        sampledItems = sorted(list(sampledItems))

        if (MAX_LINK_LENGTH<3e6):
            #If the longest link is relatively short, then no new noised link may be longer than it. This is needed to preserve the overall dataset properties (e.g. if all links >2Mbp ar deleted from the dataset)
            #for each node, I find all those that are nearby
            for node in set(sampledItems):
                pool = [el for el in items if abs(el-node)<=MAX_LINK_LENGTH] #all other ends for a link
                poolWeights = [nodeDegrees[el] for el in pool]
                poolWeights = [el/sum(poolWeights) for el in poolWeights]
                # if sum(poolWeights)!=1:
                #     raise Exception("poolWeights sum should be 1 as it is normalized")
                k = sampledItems.count(node)  # k is desired number of other ends for this node
                otherEnds = np.random.choice(pool, size=k, replace=True, p=poolWeights)
                for oe in otherEnds:
                    A = node
                    B=oe
                    RLNOISE.append(tuple([A,B]) if (A<B) else tuple([B,A]))
            iii=0
        else:
            #If longest link is long, then I divide the chomosome into 'bags' and first choose one 'bag' respecting probabilities of choosing each of them, and only then a segment from the bag
            start = min(nodeDegrees.keys())
            end = max(nodeDegrees.keys())
            bagCount = math.floor(math.log2(end-start))+1
            bags = [[] for _ in range(bagCount)]
            sortedNodes = list(sorted(nodeDegrees.keys()))
            bagSize = (end-start)/bagCount
            lastCoordinate = start+bagSize
            bagStarts = [int(start+i*bagSize) for i in range(bagCount)]
            bagInd = 0
            i=0
            while(i<len(sortedNodes)):
                node = sortedNodes[i]
                if node<=lastCoordinate:
                    bags[bagInd].append(node)
                    i+=1
                else:
                    lastCoordinate+=bagSize
                    if bagInd<len(bags)-1:
                        bagInd+=1

            
            bagTotalWeights = [(sum([nodeDegrees[node] for node in bag])/sum(nodeDegrees.values())) for bag in bags]
            bagWeights = [[nodeDegrees[node] for node in bag] for bag in bags]
            bagWeights = [[el/sum(bagWeights[i]) for el in bagWeights[i]] for i in range(len(bagWeights))]
            bagSize = int(bagSize)
            for node in sampledItems:
                minPossibleBag = max(0,math.floor((node - MAX_LINK_LENGTH)/bagSize))
                maxPossibleBag = min(bagCount-1, math.ceil((node+MAX_LINK_LENGTH)/bagSize))
                if sum(bagTotalWeights[minPossibleBag:maxPossibleBag+1])>0:
                    bagChosen = random.choices(list(range(minPossibleBag, maxPossibleBag+1)), k=1, weights=bagTotalWeights[minPossibleBag:maxPossibleBag+1])[0]
                    otherEnd = np.random.choice(bags[bagChosen], size=1, replace=True, p=bagWeights[bagChosen])[0]
                    A = node
                    B=otherEnd
                    RLNOISE.append(tuple([A,B]) if (A<B) else tuple([B,A]))
                else:
                    continue #all bags had zero weight. Should not happen often
            iii=0




        #Adding genomic effect noise
        GENOISE = []
        noiseLengths = random.choices([abs(el[1]-el[0]) for el in actualOLinks], k=gdnoiseToAdd)
        #These are the new link lengths.
        #noiseLengths is a list of link lengths that I will be creating

        allNodes = sorted(list(nodeDegrees.keys()))
        allNodesSet = set(allNodes)
        i=0
        tolerance = math.floor(math.log2(len(allNodesSet))+1) #approx. 14, nr of approximated node positions I will look in both directions to find a node that was a part of the original dataset. needed because in some approximated positions there moght be no nodes in the original dataset
        diffTol = tolerance//2
        #closeDeltas2 = [(j-diffTol)*minDiff for j in range(tolerance)]
        closeDeltas = [j*minDiff for j in range(-diffTol, tolerance-diffTol, 1)]
        iters = 0
        selfLoops=0
        while (len(GENOISE)<gdnoiseToAdd):
            #choose one end randomly
            A = random.choice(allNodes)
            if (noiseLengths[i]==0):
                i+=1
                newLink = (A, A) 
                GENOISE.append(newLink)
                continue
            sign = -1 if random.randint(0, 1) == 0 else +1
            B = A+(noiseLengths[i]*sign)
            bpBetweenAAndB = abs(A-B)
            setOfBs = set([B+delta for delta in closeDeltas if abs(delta)<bpBetweenAAndB])
            if A in setOfBs:
                raise Exception("Self loop can be created")
            existingBs = allNodesSet.intersection(setOfBs)
            if len(existingBs)>0:
                B = random.choice(list(existingBs))
                newLink = (A, B) if (A < B) else (B, A)
                GENOISE.append(newLink)
                i+=1
            iters+=1

        print(f"Created {len(GENOISE)} set with {len(set(GENOISE))} different elements in {iters} iterations, and {selfLoops=}")

        ORIGLINKS = [(A,B) for [A,B] in actualOLinks]

        #Start downsampling
        ALLLINKS = ORIGLINKS+GENOISE+RLNOISE
        DOWNSAMPLEDLINKS  = random.sample(ALLLINKS, len(ORIGLINKS))
        DOWNSAMPLEDLINKS = set(DOWNSAMPLEDLINKS) #remove duplicates
        its = 0
        while(len(DOWNSAMPLEDLINKS)<len(ORIGLINKS)):
            extraLinks =  random.sample(ALLLINKS, len(ORIGLINKS)-len(DOWNSAMPLEDLINKS))
            DOWNSAMPLEDLINKS.update(extraLinks)
            its+=1
        print(len(ORIGLINKS), "-->", len(DOWNSAMPLEDLINKS), " in ", its)

        #Now, create Universal
        for link in DOWNSAMPLEDLINKS:
            A = segmentStartToIndex[link[0]]
            B = segmentStartToIndex[link[1]]
            if A>B: A,B=B,A
            if (A,B) not in randomizedLinkBits:
                randomizedLinkBits[(A,B)]=0
            randomizedLinkBits[(A,B)]|=tisBit
        print(f"Done {tis}")
    #end loop over tissues of the chr
    R["chrValues"][ch]["links"] = [[A, B, bit] for (A, B), bit in sorted(randomizedLinkBits.items(), key=lambda x: (x[0][0], x[0][1]))]
    with open(noisyFn, 'w') as file:
        json.dump(R, file)
    print(noisyFn, "created")