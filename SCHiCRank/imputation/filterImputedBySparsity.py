#!/usr/bin/env python
"""
Filter imputed HDF5 files by sparsity (mean node degree).

For a given r (mean node degree), keeps interactions by count
such that: r = #interactions / #active_bins
Where active_bins = set of all bins participating in interactions

Usage:
    python filterImputedBySparsity.py --input-hdf5 <file> --output-hdf5 <file> --r <value> --n-bins <int>
"""

import argparse
import h5py
import numpy as np
import sys
from scipy.sparse import csr_matrix


def filter_by_mean_degree(hdf5_input, hdf5_output, r, n_bins):
    """
    Filter imputed interactions to achieve target mean node degree.

    Args:
        hdf5_input: Input HDF5 file with imputed data
        hdf5_output: Output HDF5 file for filtered data
        r: Target mean node degree (r = interactions/active_bins)
        n_bins: Number of bins in the matrix (for matrix size only)

    Returns:
        Number of interactions kept
    """
    # Read the imputed HDF5
    with h5py.File(hdf5_input, "r") as f:
        group = f["Matrix"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]

    # Create sparse matrix
    matrix_sparse = csr_matrix((data, indices, indptr), shape=(n_bins, n_bins))
    matrix_coo = matrix_sparse.tocoo()

    # Get all interactions with their counts
    interactions = np.column_stack([matrix_coo.row, matrix_coo.col, matrix_coo.data])

    # Remove diagonal (self-interactions)
    interactions = interactions[interactions[:, 0] != interactions[:, 1]]

    # Sort by count (descending)
    sorted_indices = np.argsort(interactions[:, 2])[::-1]
    interactions_sorted = interactions[sorted_indices]

    rows = interactions_sorted[:, 0].astype(int)
    cols = interactions_sorted[:, 1].astype(int)
    n_interactions = len(interactions_sorted)

    # Phase 2 prep: for each bin, find the index of its strongest interaction
    # np.minimum.at keeps the smallest index = strongest signal (sorted descending)
    bin_best = np.full(n_bins, n_interactions, dtype=int)
    all_indices = np.arange(n_interactions)
    np.minimum.at(bin_best, rows, all_indices)
    np.minimum.at(bin_best, cols, all_indices)

    # Mandatory interactions: one per bin that has any interaction (phase 2)
    has_interaction = bin_best < n_interactions
    mandatory_indices = np.unique(bin_best[has_interaction])

    # Active bins count is fixed after phase 2 (all bins with any interaction)
    n_active = int(np.sum(has_interaction))

    # Target interaction count: r * active_bins, at least n_bins/5 (phase 1)
    initial_batch = int(n_bins / 5)
    target_n = max(int(np.ceil(r * n_active)), initial_batch)
    target_n = min(target_n, n_interactions)

    # Build kept set: top target_n by strength + mandatory
    kept_mask = np.zeros(n_interactions, dtype=bool)
    kept_mask[:target_n] = True
    kept_mask[mandatory_indices] = True

    top_interactions = interactions_sorted[kept_mask]
    n_keep = int(np.sum(kept_mask))
    final_degree = n_keep / n_active if n_active > 0 else 0

    # Original active bins = same as n_active (all bins with any interaction)
    original_active_bins = n_active

    # Return stats for batch aggregation instead of printing per-cell
    stats = {
        'total_interactions': len(interactions),
        'kept_interactions': n_keep,
        'active_bins': n_active,
        'original_active_bins': original_active_bins,
        'n_bins': n_bins,
        'final_degree': final_degree,
    }

    # Create new sparse matrix with filtered interactions
    row = top_interactions[:, 0].astype(int)
    col = top_interactions[:, 1].astype(int)
    data_filtered = top_interactions[:, 2]

    # Create symmetric matrix (add both directions)
    row_sym = np.concatenate([row, col])
    col_sym = np.concatenate([col, row])
    data_sym = np.concatenate([data_filtered, data_filtered])

    # Create CSR matrix
    filtered_matrix = csr_matrix((data_sym, (row_sym, col_sym)), shape=(n_bins, n_bins))

    # Write to HDF5
    with h5py.File(hdf5_output, "w") as f:
        grp = f.create_group("Matrix")
        grp.create_dataset("data", data=filtered_matrix.data)
        grp.create_dataset("indices", data=filtered_matrix.indices)
        grp.create_dataset("indptr", data=filtered_matrix.indptr)

    return stats


def filter_batch(input_dir, output_dir, r, n_bins):
    """
    Filter all HDF5 files in input_dir, writing results to output_dir.
    Runs in a single process to avoid repeated Python/numpy/scipy startup.
    """
    import os
    import gc
    import glob

    os.makedirs(output_dir, exist_ok=True)

    hdf5_files = sorted(glob.glob(os.path.join(input_dir, "*.hdf5")))
    total = len(hdf5_files)
    if total == 0:
        print(f"  No HDF5 files found in {input_dir}")
        return 0

    filtered = 0
    failed = 0
    total_original_bins = 0
    total_kept_bins = 0
    total_original_interactions = 0
    total_kept_interactions = 0

    for i, hdf5_file in enumerate(hdf5_files):
        basename = os.path.basename(hdf5_file)
        output_hdf5 = os.path.join(output_dir, basename)

        try:
            stats = filter_by_mean_degree(hdf5_file, output_hdf5, r, n_bins)
            filtered += 1
            total_original_bins += stats['original_active_bins']
            total_kept_bins += stats['active_bins']
            total_original_interactions += stats['total_interactions']
            total_kept_interactions += stats['kept_interactions']

            if (i + 1) % 100 == 0:
                print(f"    Cell {i + 1}/{total}: "
                      f"n_bins={stats['n_bins']}, "
                      f"active_bins: {stats['original_active_bins']} -> {stats['active_bins']}, "
                      f"interactions: {stats['total_interactions']} -> {stats['kept_interactions']}")
                gc.collect()

        except Exception as e:
            failed += 1
            print(f"    Warning: Failed to filter {basename}: {e}")

    # Print summary for this chromosome
    if failed > 0:
        print(f"    Failed: {failed}")
    if filtered > 0:
        avg_orig_bins = total_original_bins / filtered
        avg_kept_bins = total_kept_bins / filtered
        avg_orig_int = total_original_interactions / filtered
        avg_kept_int = total_kept_interactions / filtered
        bin_retention = (avg_kept_bins / avg_orig_bins * 100) if avg_orig_bins > 0 else 0
        print(f"    Summary: {filtered} cells | "
              f"Avg bins: {avg_orig_bins:.0f} -> {avg_kept_bins:.0f} ({bin_retention:.0f}% retained) | "
              f"Avg interactions: {avg_orig_int:.0f} -> {avg_kept_int:.0f}")

    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Filter imputed HDF5 by mean node degree (sparsity)"
    )
    parser.add_argument("--input-hdf5", help="Input imputed HDF5 file (single-file mode)")
    parser.add_argument("--output-hdf5", help="Output filtered HDF5 file (single-file mode)")
    parser.add_argument("--input-dir", help="Input directory of HDF5 files (batch mode)")
    parser.add_argument("--output-dir", help="Output directory for filtered files (batch mode)")
    parser.add_argument("--r", type=float, required=True, help="Target mean node degree")
    parser.add_argument("--n-bins", type=int, required=True, help="Number of bins")

    args = parser.parse_args()

    try:
        if args.input_dir and args.output_dir:
            # Batch mode: process all files in directory
            n_filtered = filter_batch(
                args.input_dir, args.output_dir, args.r, args.n_bins
            )
            print(f"  ✅ Batch filtered {n_filtered} files")
            return 0
        elif args.input_hdf5 and args.output_hdf5:
            # Single-file mode (backwards compatible)
            filter_by_mean_degree(
                args.input_hdf5, args.output_hdf5, args.r, args.n_bins
            )
            return 0
        else:
            print("Error: Provide either --input-hdf5/--output-hdf5 or --input-dir/--output-dir")
            return 1

    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
