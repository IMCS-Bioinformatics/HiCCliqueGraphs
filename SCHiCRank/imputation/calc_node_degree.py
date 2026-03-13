#!/usr/bin/env python
"""
Calculate average node degree for each chromosome in a scool file.

Usage:
    python calc_node_degree.py --scool-file <file.scool>
"""

import argparse
import cooler
import h5py
import numpy as np


def list_scool_cells(scool_path):
    """List cells using h5py fallback for cooler API compatibility."""
    try:
        return cooler.fileops.list_scool_cells(scool_path)
    except (OSError, AttributeError):
        with h5py.File(scool_path, 'r') as f:
            if 'cells' in f:
                return [f'/cells/{name}' for name in f['cells'].keys()]
        raise


def calc_node_degree_per_chrom(scool_path):
    """
    Calculate average node degree for each chromosome in scool file.

    Args:
        scool_path: Path to scool file

    Returns:
        Dictionary mapping chromosome to average node degree
    """
    # Get list of all cells
    cell_list = list_scool_cells(scool_path)
    print(f"Found {len(cell_list)} cells in scool file")
    print("")

    # Get first cell to find chromosomes
    first_cell = cell_list[0]
    clr = cooler.Cooler(f"{scool_path}::{first_cell}")
    chromosomes = list(clr.chromsizes.keys())
    resolution = clr.binsize

    print(f"Resolution: {resolution}")
    print(f"Chromosomes: {len(chromosomes)}")
    print("")

    # Calculate node degree per cell per chromosome
    chrom_results = {}

    for chrom in chromosomes:
        # Store per-cell metrics for this chromosome
        cell_degrees = []
        cell_interactions = []
        cell_active_bins = []

        for cell_name in cell_list:
            clr = cooler.Cooler(f"{scool_path}::{cell_name}")

            # Get bins and pixels for this chromosome
            bins = clr.bins()[:]
            chrom_bins = bins[bins['chrom'] == chrom]
            n_bins = len(chrom_bins)

            if n_bins == 0:
                continue

            chrom_bin_ids = set(chrom_bins.index)
            pixels = clr.pixels()[:]

            # Get interactions for this chromosome (non-zero count)
            chrom_pixels = pixels[
                pixels['bin1_id'].isin(chrom_bin_ids) &
                pixels['bin2_id'].isin(chrom_bin_ids) &
                (pixels['count'] > 0)
            ]

            n_interactions = len(chrom_pixels)

            # Calculate unique bins involved in interactions
            # unique_bins = set(bin1_ids) ∪ set(bin2_ids)
            unique_bins = set(chrom_pixels['bin1_id']).union(set(chrom_pixels['bin2_id']))
            n_unique_bins = len(unique_bins)

            # Calculate degree for this cell
            if n_unique_bins > 0:
                cell_degree = n_interactions / n_unique_bins
            else:
                cell_degree = 0

            cell_degrees.append(cell_degree)
            cell_interactions.append(n_interactions)
            cell_active_bins.append(n_unique_bins)

        # Calculate statistics across cells
        if len(cell_degrees) > 0:
            mean_degree = np.mean(cell_degrees)
            std_degree = np.std(cell_degrees, ddof=1) if len(cell_degrees) > 1 else 0
            mean_interactions = np.mean(cell_interactions)
            mean_active_bins = np.mean(cell_active_bins)
        else:
            mean_degree = 0
            std_degree = 0
            mean_interactions = 0
            mean_active_bins = 0

        chrom_results[chrom] = {
            'mean_degree': mean_degree,
            'std_degree': std_degree,
            'mean_interactions': mean_interactions,
            'mean_active_bins': mean_active_bins,
            'n_cells': len(cell_degrees)
        }

    return chrom_results, len(cell_list)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate average node degree per chromosome in scool file"
    )
    parser.add_argument("--scool-file", required=True, help="Input scool file")

    args = parser.parse_args()

    try:
        results, n_cells = calc_node_degree_per_chrom(args.scool_file)

        print("=" * 80)
        print(f"Node Degree per Chromosome ({n_cells} cells)")
        print("=" * 80)
        print(f"{'Chromosome':<15} {'Active Bins':<12} {'Interactions':<15} {'Mean Degree':<15} {'StdDev':<10}")
        print("-" * 80)

        total_bins = 0
        total_interactions = 0
        all_degrees = []

        for chrom in sorted(results.keys(), key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').zfill(2))):
            data = results[chrom]
            if data['n_cells'] > 0:
                print(f"{chrom:<15} {data['mean_active_bins']:<12.1f} {data['mean_interactions']:<15.1f} {data['mean_degree']:<15.3f} {data['std_degree']:<10.3f}")
                total_bins += data['mean_active_bins']
                total_interactions += data['mean_interactions']
                all_degrees.append(data['mean_degree'])

        print("-" * 80)
        overall_mean = np.mean(all_degrees) if len(all_degrees) > 0 else 0
        overall_std = np.std(all_degrees, ddof=1) if len(all_degrees) > 1 else 0
        print(f"{'OVERALL':<15} {total_bins:<12.1f} {total_interactions:<15.1f} {overall_mean:<15.3f} {overall_std:<10.3f}")
        print("=" * 80)
        print(f"\nNote: Node degree per cell = #interactions / #active_bins")
        print(f"      Active bins = set of all bins participating in interactions")
        print(f"      Statistics: mean ± stdev across {n_cells} cells per chromosome")

        return 0

    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
