# Gene ontology enrichment analysis
This directory has scripts that are used to perform GOEA

runGOEA.py generates files with results. Note, it includes hard-coded paths to data files. So, it may be necessary to execute the script from the root directory of this project.

Directory ontologyData contains 2 files that allow to map human genes to chromosome segments and GO terms to genes.

### What is GOEA
*This is from our paper*
To analyze our results, we performed Gene Ontology Enrichment Analysis (GOEA) using the GOATools library. This analysis seeks to identify Gene Ontology terms that are overrepresented in a study gene list compared to a population gene list.

Initially, we assigned genes to nodes in the graph G. Nodes located on chromosomes where a particular gene is expressed are assigned that gene *hg19_genes_w-go.txt was used*. Subsequently, we calculated topological elements, such as C3 and S(k), identified a subset of nodes that form these topological elements, and collated a set of genes present in at least one node of the topological elements. We then conducted the GOEA. *Class Biolog from GOEAa.py module*

We executed several analyses, comparing different topological elements. For instance, we collected all C3 elements with a degree of 2 or more and took the genes present in at least one such C3 to form the study gene list. This list was then compared to the population gene list, which consisted of all genes found in at least one link in the graph.

In addition, we compared genes of supports S(k) to genes of C3 and to genes of all links. We enumerated all supports S(k), identified a list of genes present in at least one node of any support S(k) to form the study gene list, and compared it with the list of genes found in at least one node of any C3, which served as the population gene list. Different values of support degrees k and C3 degrees were used in different runs to detect any Gene Ontology terms that appear more frequently in supports than in C3. Genes from all graph nodes were also used to create the population gene list and run the analysis, comparing S(k) genes to them.

This strategy enabled us to identify Gene Ontology terms that are disproportionately represented in nodes of supports compared to all nodes of the graph.

