# SCHiCRank Pipeline

Full pipeline for single-cell Hi-C chromatin motif ordering.
Entry point: `run_full_100kb_pipeline.sh`

## Input data

```
sourceData/
    nagano_10kb_cell_types.scool          # Single-cell Hi-C at 10kb resolution (scool/HDF5)
    nagano_assoziated_cell_types.txt      # Cell name -> cell cycle phase mapping (TSV)
    chrom_sizes.txt                       # Chromosome -> length in bp (TSV)
```

## Environment setup

```
setup_schicluster_env.sh
    conda create schicluster_env4 (python=3.6.8)
    pip install scHiCluster (from github:zhoujt1994/scHiCluster)
    pip install cooler, h5py, networkx, kneed, scipy, matplotlib, opencv-python
```

## Pipeline call tree

```
run_full_100kb_pipeline.sh
│
├── STAGE 1: Resolution change (10kb -> 100kb)
│   │
│   └── resolutionChange/coarsen_scool.py --input-scool ... --output-scool ... --factor 10
│           main()
│               coarsen_scool(input_scool, output_scool, factor)
│                   cooler.Cooler()                        # load each cell from scool
│                   cooler.coarsen_cooler()                # bin merging by factor k
│                   cooler.create_scool()                  # write merged cells to new scool
│       Output: resolutionChange/coarsened_k10_*.scool
│
├── STAGE 2: Imputation
│   │
│   └── imputation/run_full_imputation.sh <coarsened.scool> ALL <r>
│       │
│       ├── 2.1  exportScoolToText.py --scool ... --chromosomes ALL
│       │            main()
│       │                generate_chrom_sizes(scool, output_file)
│       │                    cooler.Cooler()                # read bin table for chrom lengths
│       │                export_cells_for_chromosome(scool, chrom, output_dir)
│       │                    list_scool_cells(scool)       # h5py: list /cells/ group
│       │                    cooler.Cooler()               # fetch pixels per cell per chrom
│       │                    write tab-separated text      # bin1 bin2 count
│       │        Output: imputation/raw/{chrom}/{cell}.txt
│       │                imputation/chrom_sizes.txt
│       │
│       ├── 2.2  hicluster impute-cell (external, per cell × per chrom)
│       │            --pad 1 --std 1 --rp 0.5
│       │            Random walk with restart imputation (scHiCluster)
│       │        Output: imputation/imputed/{chrom}/{cell}_{chrom}_pad1_std1_rp0.5_sqrtvc.hdf5
│       │
│       ├── 2.3  filterImputedBySparsity.py --input-dir ... --output-dir ... --r <r>
│       │            main()
│       │                filter_batch(input_dir, output_dir, r, n_bins)
│       │                    filter_by_mean_degree(hdf5_in, hdf5_out, r, n_bins)
│       │                        h5py: load dense matrix
│       │                        keep top N = ceil(r * active_bins) interactions
│       │                        save as sparse CSR in HDF5
│       │        Output: imputation/filtered_r{R}/{chrom}/*.hdf5
│       │
│       ├── 2.4  convertImputedToScool.py --input-dir ... --output ...
│       │            main()
│       │                create_bins_dataframe(chrom_sizes, resolution)
│       │                process_imputed_cells(imputed_dir, bins_df, chromosomes)
│       │                    read_imputed_hdf5(hdf5_file)          # load sparse CSR
│       │                    map bin indices with chromosome offsets
│       │                create_scool_from_cells(cell_coolers, bins_df, output_scool)
│       │                    cooler.create_cooler()        # per cell
│       │                    cooler.create_scool()         # combine into multi-cell
│       │        Output: imputation/imputed_pad1_std1_rp0.5_r{R}_coarsened_k10_*.scool
│       │
│       └── 2.5  calc_node_degree.py --scool ... (verification)
│                    main()
│                        calc_node_degree_per_chrom(scool)
│                            list_scool_cells(scool)
│                            cooler.Cooler()               # compute mean degree per chrom
│                    Prints: table of mean node degrees
│
├── STAGE 3: Motif ordering
│   │
│   ├── 3.0  generateCellPhaseInfo.py --scool-file ... --cell-type-file ... --work-dir ...
│   │            main()
│   │                generate_cell_phase_info(scool, cell_type_file, output, work_dir)
│   │                    MulticoolProcessor(scool)
│   │                        MulticoolProcessor.__init__()     # h5py: list /cells/
│   │                        MulticoolProcessor.getCellNames() # extract cell names
│   │                    match cell names to phase annotations
│   │                    pickle.dump({cell_names, cell_phase})
│   │        Output: {work_dir}/cellAndPhaseInfo.pkl
│   │
│   ├── 3.1  runProcessCoolDataset.py -k 1 --chromosomes ALL --fn <scool> --work-dir ...
│   │            main()
│   │                processOriginalCoolDataset.process_cells(k, chromosomes, fn, ...)
│   │                    MulticoolProcessor(fn)
│   │                        MulticoolProcessor.readCell(name)
│   │                            CoolProcessor(path, cellName)
│   │                                CoolProcessor.reduceResolution(k)   # if k > 1
│   │                    for each chromosome:
│   │                        CoolProcessor.getAllInteractionsWithLoci(ch)
│   │                            CoolProcessor.__setAllInteractionsWithLoci(ch)
│   │                                CoolProcessor.getInteractions(ch)
│   │                                    CoolProcessor.__setInteractions(ch)
│   │                                translate bin indices -> genomic loci
│   │                        combineDicts(cell_links, link_cells)
│   │                        pickle.dump({cell_links, link_cells, ...})
│   │        Output: {work_dir}/{postfix}-{chrom}-{resolution}.pkl  (×20 chromosomes)
│   │
│   ├── 3.2  runCreateCliques.py --chromosomes ALL --postfix ... --work-dir ...
│   │            main()
│   │                for each chromosome:
│   │                    createCliqueDatafiles.createCliquePickles(baseFn, resFn)
│   │                        pickle.load(interaction pkl)
│   │                        for each cell:
│   │                            networkx.Graph() from cell_links
│   │                            networkx.find_cliques()           # all maximal cliques
│   │                            filter by size: K3, K4, K5, K6, K7, K8
│   │                        build reverse mapping: clique -> set(cells)
│   │                        pickle.dump({cell_cliques, clique_cells})
│   │        Output: {work_dir}/{postfix}-{chrom}-{resolution}-cliques.pkl  (×20)
│   │
│   ├── 3.3  runCreatePairwiseSimilarities.py --chromosomes ALL --work-dir ...
│   │            main()
│   │                for each chromosome:
│   │                    createPairwiseSimilarities.callPairwiseSimilarites(clique_pkl)
│   │                        pickle.load(clique pkl)
│   │                        for each motif type (K3..K8) × length (alllengths, long):
│   │                            for each clique shared by ≥2 cells:
│   │                                if long: filter span ≥ 2Mb
│   │                                count all pairwise cell combinations
│   │                            sort pairs by frequency descending
│   │                            write CSV (Item1, Item2, Frequency, phases, ...)
│   │        Output: {work_dir}/pairwiseSimilarities/{K}-{length}/*.csv  (×20 per motif)
│   │
│   ├── 3.4  [Optional] createInteractionPairwiseSimilarities.py --work-dir ... --config ...
│   │            main()
│   │                generate_csvs(work_dir, config)
│   │                    for each chromosome:
│   │                        pickle.load(interaction pkl)
│   │                        build sparse incidence matrix (interactions × cells)
│   │                        co-occurrence = A.T @ A
│   │                        keep top 20 partners per cell
│   │                        write CSV for alllengths and long (span ≥ 2Mb)
│   │        Output: {work_dir}/pairwiseSimilarities/interaction-{length}/*.csv
│   │
│   ├── 3.5  [Optional] PageRank ordering (for each motif config):
│   │   │
│   │   ├── runPageRankFilter.py --input-dir <pairwiseSimilarities/K3-alllengths/> --k 5
│   │   │       main()
│   │   │           runSCHiCRank.run_pagerank_filter(input_dir, label, ...)
│   │   │               build_full_neighbor_map(input_dir)
│   │   │                   for each chrom CSV:
│   │   │                       pd.read_csv() -> stack both directions
│   │   │                       sort by (cell, -frequency)
│   │   │                       store: {cell: [(partner, freq), ...]}
│   │   │                   pickle.dump(full_neighbor_map)
│   │   │               while active_cells > threshold:
│   │   │                   for each chromosome:
│   │   │                       build nx.DiGraph: cell -> top K active neighbors
│   │   │                       nx.pagerank(G)
│   │   │                       collect per-cell scores
│   │   │                   trimmed aggregation: sort scores, drop 2 best + 2 worst, sum
│   │   │                   sort cells by aggregate score descending
│   │   │                   KneeLocator(curve=convex, direction=decreasing)
│   │   │                   remove cells above elbow -> record (cell, iteration, score)
│   │   │               save remaining cells
│   │   │       Output: {work_dir}/final_active_cells_{label}.csv
│   │   │
│   │   └── runCalculateTau.py --input-csv final_active_cells_*.csv --motif-config ...
│   │           main()
│   │               pd.read_csv(input_csv)
│   │               for each permutation of phases:
│   │                   calculate_tau_for_phase_ordering(df, phase_order)
│   │                       assign ideal ranks from phase order
│   │                       scipy.stats.kendalltau(actual_ranks, ideal_ranks)
│   │               report best tau + phase ordering
│   │           Output: append to {work_dir}/tau_scores.txt
│   │
│   ├── 3.6  Motif-count ordering (for each K3..K8 × alllengths/long):
│   │   │
│   │   └── calculateTauFromMotifCounts.py --motif-type K3-alllengths --work-dir ...
│   │           main()
│   │               load_cell_phases(cellAndPhaseInfo.pkl)
│   │               collect_motif_counts_per_chromosome(chromosomes, postfix, resolution)
│   │                   for each chrom: pickle.load(clique pkl) -> count motifs per cell
│   │               rank_cells_by_motif_count(counts, motif_type, min_count=10)
│   │                   per-chrom rank via argsort -> median rank across chroms
│   │               calculate_best_tau(ranked_cells, cell_phases)
│   │                   for each of 120 phase permutations:
│   │                       calculate_tau_for_phase_ordering()
│   │                           scipy.stats.kendalltau()
│   │           Output: append to {work_dir}/tau_scores.txt
│   │
│   └── 3.7  plotTauScores.py --input tau_scores.txt --output tau_comparison.png
│               main()
│                   parse_tau_scores(log_file)        # regex parse tau_scores.txt
│                   plot_tau_scores(results, output)   # matplotlib bar chart
│           Output: {work_dir}/tau_comparison.png
│
└── Done. Results in motifOrdering/{postfix}/
```

## Output summary

```
motifOrdering/{postfix}/
    cellAndPhaseInfo.pkl                                  # cell metadata
    {postfix}-{chrom}-{resolution}.pkl                    # interactions (×20)
    {postfix}-{chrom}-{resolution}-cliques.pkl            # cliques K3-K8 (×20)
    pairwiseSimilarities/{K}-{length}/*.csv               # cell-pair co-occurrence
    final_active_cells_*.csv                              # PageRank orderings
    tau_scores.txt                                        # all Kendall's tau results
    tau_comparison.png                                    # visualization
```

## File inventory

```
run_full_100kb_pipeline.sh              # main orchestrator
setup_schicluster_env.sh                # conda environment setup

sourceData/
    nagano_10kb_cell_types.scool        # input data (not included, 741MB)
    nagano_assoziated_cell_types.txt    # cell phase annotations
    chrom_sizes.txt                     # chromosome sizes

resolutionChange/
    coarsen_scool.py                    # 10kb -> 100kb bin merging

imputation/
    run_full_imputation.sh              # imputation orchestrator
    exportScoolToText.py                # scool -> text for scHiCluster
    filterImputedBySparsity.py          # sparsify imputed matrices to target degree
    convertImputedToScool.py            # filtered HDF5 -> scool
    calc_node_degree.py                 # verify mean node degree

motifOrdering/
    MulticoolProcessor.py              # scool reader class
    CoolProcessor.py                   # single-cell Hi-C matrix class
    generateCellPhaseInfo.py           # extract cell phase metadata
    runProcessCoolDataset.py           # CLI: scool -> per-chrom interaction pkls
    processOriginalCoolDataset.py      # core: extract interactions per cell per chrom
    runCreateCliques.py                # CLI: find K3-K8 cliques
    createCliqueDatafiles.py           # core: networkx.find_cliques per cell
    runCreatePairwiseSimilarities.py   # CLI: clique co-occurrence CSVs
    createPairwiseSimilarities.py      # core: count shared cliques between cell pairs
    createInteractionPairwiseSimilarities.py  # co-occurrence from raw interactions
    runPageRankFilter.py               # CLI: iterative PageRank cell ordering
    runSCHiCRank.py                    # core: k-NN graph + PageRank + elbow detection
    runCalculateTau.py                 # evaluate PageRank ordering via Kendall's tau
    calculateTauFromMotifCounts.py     # evaluate motif-count ordering via Kendall's tau
    plotTauScores.py                   # bar chart of tau scores
```
