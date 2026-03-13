# SCHiCRank

A pipeline for ordering single-cell Hi-C (scHi-C) cells by cell-cycle pseudotime using topological motifs in chromatin contact graphs.

## What it does

SCHiCRank takes a dataset of single-cell Hi-C contact matrices (in `.scool` format) and produces a pseudotime ordering of cells that reflects cell-cycle progression. The approach works by:

1. **Resolution coarsening** -- Merging genomic bins (e.g. 10kb to 100kb) to increase contact density per cell.
2. **Imputation** (optional) -- Applying random-walk-with-restart smoothing (via scHiCluster) to fill in missing contacts, then filtering to a controlled sparsity level by mean node degree. This is a standard step in scHi-C data processing.
3. **Motif detection** -- For each cell and each chromosome, building a graph where nodes are genomic bins and edges are detected contacts, then enumerating all maximal cliques of sizes 3--8 (K3--K8). These cliques represent local topological motifs in the chromatin structure.
4. **Pairwise similarity** -- For each pair of cells, counting how many motif instances they share across all chromosomes.
5. **Cell ordering** -- Two methods are provided:
   - **SCHiCRank (PageRank-based):** For each chromosome, a directed k-nearest-neighbor graph is built from motif co-occurrence. PageRank centrality is computed per chromosome and aggregated across chromosomes via a trimmed sum. An elbow is detected in the sorted centrality curve, and cells above the elbow are removed and placed into the ordering. This repeats iteratively, peeling successive clusters.
   - **Motif-count ordering:** Cells are ranked by their total motif count per chromosome, then aggregated across chromosomes by median rank.
6. **Evaluation** -- The resulting ordering is compared against known cell-cycle phase labels using Kendall's tau, maximized over all cyclic permutations of the five phases (G1, early-S, late-S/G2, pre-M, post-M).

## Quick start

### Prerequisites

- Linux (imputation requires scHiCluster, which only runs on Linux)
- Conda

### Setup

```bash
# Create environment and install scHiCluster + dependencies
cd SCHiCRank
./setup_schicluster_env.sh
```

### Input data

Download from [Zenodo (10.5281/zenodo.4308298)](https://zenodo.org/records/4308298):
- `nagano_10kb_cell_types.scool` (741 MB) -- place in `sourceData/`. Note that for demo purposes this repo contains a version of the scool file with just 10 cells. Please manually upload the full file and replace the provided file.
- Cell type annotations are already included in `sourceData/nagano_assoziated_cell_types.txt`

### Run

```bash
cd SCHiCRank
./run_full_100kb_pipeline.sh
```

This runs the full pipeline: coarsen 10kb to 100kb, impute with r=3.0, detect cliques, compute pairwise similarities, and evaluate orderings. Results are written to `motifOrdering/imputed100k_r3.0/`.

See [pipeline.md](pipeline.md) for a detailed call tree of every script and function.

## Repository structure

```
run_full_100kb_pipeline.sh         # Main pipeline orchestrator
setup_schicluster_env.sh           # Conda environment setup

sourceData/                        # Input data
    nagano_10kb_cell_types.scool   # scHi-C dataset (not included, download above)
    nagano_assoziated_cell_types.txt
    chrom_sizes.txt

resolutionChange/                  # Stage 1: bin coarsening
    coarsen_scool.py

imputation/                        # Stage 2: scHiCluster imputation + sparsity filtering
    run_full_imputation.sh
    exportScoolToText.py
    filterImputedBySparsity.py
    convertImputedToScool.py
    calc_node_degree.py

motifOrdering/                     # Stage 3: motif analysis + cell ordering
    CoolProcessor.py               # Single-cell contact matrix class
    MulticoolProcessor.py           # Multi-cell scool reader
    generateCellPhaseInfo.py        # Extract cell phase metadata
    runProcessCoolDataset.py        # scool -> per-chromosome interaction pickles
    processOriginalCoolDataset.py
    runCreateCliques.py             # Find K3-K8 cliques per cell
    createCliqueDatafiles.py
    runCreatePairwiseSimilarities.py  # Count shared cliques between cell pairs
    createPairwiseSimilarities.py
    runPageRankFilter.py            # Iterative PageRank cell ordering (SCHiCRank)
    runSCHiCRank.py
    runCalculateTau.py              # Evaluate ordering quality
    calculateTauFromMotifCounts.py  # Motif-count-based ordering + evaluation
    plotTauScores.py                # Visualization
```

## Reference dataset

The example dataset is from Nagano et al. (2017), containing 1,078 mouse embryonic stem cells with known cell-cycle phases, provided as a `.scool` file by Wolff et al. ([Zenodo](https://zenodo.org/records/4308298)).


