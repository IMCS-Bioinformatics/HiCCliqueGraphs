# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 19:00:00 2023

@author: asizo
"""
from universal import *
from topologicalFeatures import *
import json
import csv
from pathlib import Path

from intervaltree import Interval, IntervalTree

from goatools.utils import read_geneset
from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.go_enrichment import GOEnrichmentStudy


        
class GOEATool:
    def __init__(self, fin_obo, hg19_genes_w_go, dx=0, annoDir = "annos", pv=0.05):
        # fin_obo - file name for go-basic.obo file - DAG containing HPO terms
        # hg19_genes_w-go - hg19_genes_w-go.txt file with path. It has all genes, their terms and locations for all chromosomes. necessarty for Gene to segments mapping
        #this class should be created only once and then the method GOEA should be called for each analysis
        # import collections as cx
        # from goatools.utils import read_geneset
        # from goatools.obo_parser import GODag
        # from goatools.anno.idtogos_reader import IdToGosReader
        # from goatools.go_enrichment import GOEnrichmentStudy
        # these imports are required
        # https://github.com/tanghaibao/goatools.git
        self.godag = GODag(fin_obo)
        self.pv = pv
        self.dx = dx
        
        self.curResultAll = None #all terms and their p-values in the most recent analysis
        self.curResultSig = None #significant GO terms in the most recent GOE analysis
        self.annoDir = annoDir
        
        self.makeIntervalTreesForAllChrs(hg19_genes_w_go, self.dx) #make chromosome interval trees with genes, sets self.intTrees

        self.curAccumulatedResults = []
        
        
    
    def makeIntervalTreesForAllChrs(self, hg19_genes_w_go, dx=0):
        # makes interval trees for all chromosomes to be able to map genes to segments. 
        # also creates gene annotation files for goatools library
        chrs = ["chr"+str(i) for i in range(1,23,1)]
        chrs.append("chrX")
        chrs = set(chrs)
    
        ttt = {ch: IntervalTree() for ch in chrs} #dict of interval trees
        
        geneTerms = {ch: dict() for ch in chrs}
    
        with open(hg19_genes_w_go, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                header = row
                # header == Gene stable ID,Chromosome/scaffold name,Gene start (bp),Gene end (bp),Gene name,GO term accession,GO term name
                break
            for row in reader:
                # row == e.g., ENSG00000196431,22,27017928,27026636,CRYBA4,GO:0043010,camera-type eye development
                curChr = "chr"+row[1]
                if curChr not in chrs: continue
                #this row is about a chromosome that is a real chromosome, not outside of chromosomes.
                #current interval is from [ a:=row[2] to b:=row[3]+1 ) 
                # a and b are str at first, transforming them to int
                a = int(row[2])
                b = int(row[3])+1 # +1 because intervals in IntervalTree are including start and excluding end
                (a,b) = sorted((a,b))
                # dx should be added to both in order to extend the region by dx to every side
                a-=dx
                b+=dx
                gene = row[4]
                ttt[curChr][a:b] = gene #to this segment I add the gene!!!
                
                #now add go term to the gene to make annotation dictionary
                goTerm = row[5]
                if gene not in geneTerms[curChr]: geneTerms[curChr][gene] = set()
                geneTerms[curChr][gene].add(goTerm)
        
        #now make annotation folder with annotation files for every chromosome and for allChrs together
        p = Path('annos')
        p.mkdir(parents=True, exist_ok=True)
        allChrs = []
        for ch in chrs:
            nfn = "annos/{}-annos.tsv".format(ch)
            L = []
            for gene in geneTerms[ch].keys():
                strTerms = ';'.join(geneTerms[ch][gene])
                row = [gene, strTerms]
                L.append(row)
            with open(nfn, 'w', newline='') as file:
                writer = csv.writer(file, delimiter = '\t')
                writer.writerows(L)
            allChrs.extend(L)
        #write one file for all chromosomes together
        nfn = "annos/chrAll-annos.tsv"
        with open(nfn, 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerows(allChrs)
        self.intTrees = ttt        
        return ttt
    
    def getListOfGenesFromListOfSegments(self, listOfSegments, ch, all=""):
        #getsListOfGenesFromListOfSegments. self.intTrees must be set
        chrGenes = set()
        #if all=="all", then adding all possible this chromosome's genes. No matter if there's a segment in our graph or not
        if all=="all":
            genes = list(el.data for el in self.intTrees[ch][0:500e6])
            chrGenes.update(genes)
        else:
            #normal case, chromosome's filtered segmetns genes are searched
            for seg in listOfSegments:
                genes = list( el.data for el in self.intTrees[ch][int(seg[0]):int(seg[1])])
                chrGenes.update(genes)

        return sorted(list(chrGenes))    
    
            
    def GOEA(self, stuObject, popObject=None, ch=None, what="something"):
        # sets self.curResultAll and curResultSig
        # returns count of significant "enriched" GO terms
        #accepts 2 objects as parameters - one with population graph and the other with study graph
        # objects must have getListOfSegmentLociSegments() method implemented
        # population object can be None, then the object that was used the previous time will be used again
        #what is a text description of what is being compared
        #set popObject to "all" to compare to all chromsome's genes
        
        if popObject is None:
            #try to use the same object as the previous time
            try:
                popObject = self.prevPopObject
                listOfPopulationGenes = self.prevListOfPopulationGenes
                if ch is None: ch = stuObject.ch                
                #and calculate study genes the normal way, from the very beginning
                listOfStudySegments = stuObject.getListOfSegmentLociSegments()
                listOfStudyGenes = self.getListOfGenesFromListOfSegments(listOfStudySegments, ch)
            except AttributeError:
                raise ValueError("No previous popObject or chromosome found. Call GOEATool.GOEA method with all parameters")
        else:       
            #calculate both sets of genes again
            try:
                if ch is None: ch = stuObject.ch #same chromsome as study object

                if popObject=="all": #compare to all possible genes in this chromosome
                    listOfPopulationGenes = self.getListOfGenesFromListOfSegments([], ch, "all")
                else:
                    listOfPopulationSegments = popObject.getListOfSegmentLociSegments()
                    listOfPopulationGenes = self.getListOfGenesFromListOfSegments(listOfPopulationSegments, ch)
                    
                listOfStudySegments = stuObject.getListOfSegmentLociSegments()
                listOfStudyGenes = self.getListOfGenesFromListOfSegments(listOfStudySegments, ch)

                
            except AttributeError:
                raise ValueError("object passed to GOEA method of GOEATool class that has no getListOfSegmentLociSegments() method")
            

        self.prevListOfPopulationGenes = listOfPopulationGenes
        self.prevCh = ch
        self.prevPopObject = popObject
        #at this point listOfStudyGenes and listOfPopulationGenes are set. 
        
        #write gene lists to file
        with open("pop-tmp.tsv", 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            for row in listOfPopulationGenes:
                writer.writerow([row]) 
        with open("stu-tmp.tsv", 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            for row in listOfStudyGenes:
                writer.writerow([row]) 
        
        study_ids = read_geneset("stu-tmp.tsv")
        population_ids = read_geneset("pop-tmp.tsv")
        fin_anno = self.annoDir+"/"+ch+"-annos.tsv"
        
        annoobj = IdToGosReader(fin_anno, godag=self.godag)
        id2gos = annoobj.get_id2gos()
        goeaobj = GOEnrichmentStudy(
            population_ids,
            annoobj.get_id2gos(),
            self.godag,
            methods=['fdr_bh'],
            pvalcalc='fisher_scipy_stats',
            prt=None,)
        
        self.curResultAll = goeaobj.run_study(study_ids, prt=None)
        self.curResultSig = [r for r in self.curResultAll if r.p_fdr_bh < self.pv]
        
        E=sum(1 for r in self.curResultSig if r.enrichment=='e') #count of enriched GO terms
        self.curResultSig = [r for r in self.curResultSig if r.enrichment =='e']
        
        #add header to each row with DS, what and chr
        if E>0: 
            curAppendable = [[stuObject.DS, what, ch]+str(row).split('\t') for row in self.curResultSig]
            self.curAccumulatedResults.extend([row for row in curAppendable]) #add results to accumulated list
        
        #goeaobj.wr_xlsx("nbt3102_geneids.xlsx", self.curResultSig)
        
        return E
    
    def reset(self):
        #clears curAccumulatedResults
        self.curAccumulatedResults = []
    
    def dump(self, fn="mostRecentlyDumpedGOEA.csv"):
        #writes curAccumulatedResults to a file
        header = ["GO", "NS", "enrichment", "name", "ratio_in_study", "ratio_in_pop", "p_uncorrected", "depth", "study_count", "p_fdr_bh", "study_items"]
        header = ["DS", "What", "chr"] + header
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(header)
            for row in self.curAccumulatedResults:
                writer.writerow(row)
        return
    

