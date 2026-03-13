"""
Generate pairwiseSimilarities CSVs from interaction pkl files.

For each chromosome pkl, computes cell-cell co-occurrence via sparse A.T @ A,
then writes CSVs in the same format as createPairwiseSimilarities.py.
Produces two variants: interaction-alllengths and interaction-long (>=2Mb).

Usage:
    python createInteractionPairwiseSimilarities.py --work-dir motifOrdering/raw100k --config raw100k
"""

import argparse
import os
import csv
import pickle
import numpy as np
from scipy.sparse import csr_matrix


ALL_CHROMOSOMES = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chrX",
]


def generate_csvs(work_dir, config):
    for chrom in ALL_CHROMOSOMES:
        pkl_path = os.path.join(work_dir, f"{config}-{chrom}-100000.pkl")
        if not os.path.exists(pkl_path):
            print(f"  Skipping {chrom} — no pkl")
            continue

        print(f"  Loading {chrom}...", end="", flush=True)
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
        print(" done.", flush=True)

        n_cells = len(data["cell_IDs"])
        index_to_name = data["index_to_name"]
        index_to_type = data["index_to_type"]
        resolution = data["resolution"]

        # Build sparse incidence matrices
        rows_all, cols_all = [], []
        rows_long, cols_long = [], []
        inter_all_idx = 0
        inter_long_idx = 0

        for (b1, b2), cells in data["link_cells"].items():
            if len(cells) < 2:
                continue
            rows_all.extend([inter_all_idx] * len(cells))
            cols_all.extend(cells)
            inter_all_idx += 1
            if b2 - b1 >= 2_000_000:
                rows_long.extend([inter_long_idx] * len(cells))
                cols_long.extend(cells)
                inter_long_idx += 1

        for motif_length, rows, cols, n_inter in [
            ("alllengths", rows_all, cols_all, inter_all_idx),
            ("long", rows_long, cols_long, inter_long_idx),
        ]:
            result_dir = os.path.join(work_dir, f"pairwiseSimilarities/interaction-{motif_length}/")
            os.makedirs(result_dir, exist_ok=True)
            result_fn = os.path.join(result_dir,
                f"pairwiseSimilarities-{config}-{chrom}-{resolution}-interaction-{motif_length}.csv")

            if os.path.exists(result_fn):
                print(f"    {motif_length}: already exists, skipping")
                continue

            if n_inter < 1 or len(rows) == 0:
                # Write empty CSV with header
                with open(result_fn, "w", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=[
                        "Item 1", "Item 2", "Frequency",
                        "cell1_name", "cell2_name", "cell1_phase", "cell2_phase", "same?", "type"
                    ])
                    writer.writeheader()
                print(f"    {motif_length}: no interactions, wrote empty CSV")
                continue

            # Compute co-occurrence matrix
            ones = np.ones(len(rows), dtype=np.int32)
            A = csr_matrix((ones, (rows, cols)), shape=(n_inter, n_cells))
            cooc = (A.T @ A).toarray()
            np.fill_diagonal(cooc, 0)

            # For each cell, keep top TOP_K partners by co-occurrence
            TOP_K = 20
            keep_pairs = set()
            for cell_id in range(n_cells):
                row = cooc[cell_id]
                nonzero_count = np.count_nonzero(row)
                if nonzero_count == 0:
                    continue
                n_take = min(TOP_K, nonzero_count)
                top_partners = np.argpartition(-row, n_take)[:n_take]
                for partner in top_partners:
                    a, b = min(cell_id, int(partner)), max(cell_id, int(partner))
                    keep_pairs.add((a, b))

            if not keep_pairs:
                with open(result_fn, "w", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=[
                        "Item 1", "Item 2", "Frequency",
                        "cell1_name", "cell2_name", "cell1_phase", "cell2_phase", "same?", "type"
                    ])
                    writer.writeheader()
                print(f"    {motif_length}: no pairs after top-{TOP_K} filter, wrote empty CSV")
                continue

            r_idx = np.array([p[0] for p in keep_pairs])
            c_idx = np.array([p[1] for p in keep_pairs])
            freqs = cooc[r_idx, c_idx]

            # Sort by frequency descending
            order = np.argsort(-freqs)
            r_idx = r_idx[order]
            c_idx = c_idx[order]
            freqs = freqs[order]

            with open(result_fn, "w", newline="") as csvfile:
                fieldnames = ["Item 1", "Item 2", "Frequency",
                              "cell1_name", "cell2_name", "cell1_phase", "cell2_phase", "same?", "type"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for i in range(len(r_idx)):
                    c1, c2 = int(r_idx[i]), int(c_idx[i])
                    phase1 = index_to_type[c1]
                    phase2 = index_to_type[c2]
                    writer.writerow({
                        "Item 1": c1,
                        "Item 2": c2,
                        "Frequency": int(freqs[i]),
                        "cell1_name": index_to_name[c1],
                        "cell2_name": index_to_name[c2],
                        "cell1_phase": phase1,
                        "cell2_phase": phase2,
                        "same?": phase1 == phase2,
                        "type": f"{phase1}+{phase2}",
                    })

            print(f"    {motif_length}: {len(r_idx)} pairs written")


def main():
    parser = argparse.ArgumentParser(
        description="Generate pairwiseSimilarities CSVs from interaction pkl files."
    )
    parser.add_argument("--work-dir", type=str, required=True,
                        help="Config work directory (e.g., motifOrdering/raw100k)")
    parser.add_argument("--config", type=str, required=True,
                        help="Config name (e.g., raw100k)")
    args = parser.parse_args()

    print(f"Generating interaction pairwiseSimilarities for {args.config}...")
    generate_csvs(args.work_dir, args.config)
    print("Done.")


if __name__ == "__main__":
    main()
