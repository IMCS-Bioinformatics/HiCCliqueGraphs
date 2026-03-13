#!/usr/bin/env python
"""
Export cells from scool file to text format for scHiCluster imputation.

Usage:
    python exportScoolToText.py --scool-file <path> --chrom <chr> --output-dir <dir>
"""

import argparse
import cooler
import h5py
import os
import sys


def list_scool_cells(scool_path):
    """List cells using h5py fallback for cooler API compatibility."""
    try:
        return cooler.fileops.list_scool_cells(scool_path)
    except (OSError, AttributeError):
        with h5py.File(scool_path, 'r') as f:
            if 'cells' in f:
                return [f'/cells/{name}' for name in f['cells'].keys()]
        raise


def export_cells_for_chromosome(scool_path, chrom, output_dir, max_cells=None):
    """
    Export all cells from scool file for a specific chromosome to text format.

    Args:
        scool_path: Path to input scool file
        chrom: Chromosome name (e.g., 'chr1')
        output_dir: Directory to write text files
        max_cells: Maximum number of cells to export (None = all cells)

    Returns:
        Number of cells exported
    """
    os.makedirs(output_dir, exist_ok=True)

    # Get list of all cells
    cell_list = list_scool_cells(scool_path)
    print(f"Found {len(cell_list)} cells in scool file")

    # Limit number of cells if specified
    if max_cells is not None and max_cells < len(cell_list):
        cell_list = cell_list[:max_cells]
        print(f"Limiting to first {max_cells} cells for testing")

    cells_exported = 0

    # Get resolution from first cell
    first_cell = cell_list[0]
    clr_temp = cooler.Cooler(f"{scool_path}::{first_cell}")
    resolution = clr_temp.binsize

    for i, cell_name in enumerate(cell_list):
        if i % 10 == 0:
            print(f"Exporting cell {i+1}/{len(cell_list)}: {cell_name}")

        # Open the cell (cell_name already includes /cells/ prefix)
        clr = cooler.Cooler(f"{scool_path}::{cell_name}")

        # Get pixels and bins
        pixels = clr.pixels()[:]
        bins = clr.bins()[:]

        # Filter for the specific chromosome
        chrom_bins = bins[bins['chrom'] == chrom].copy()
        if len(chrom_bins) == 0:
            continue

        chrom_bin_ids = set(chrom_bins.index)

        # Filter pixels for this chromosome
        chrom_pixels = pixels[
            pixels['bin1_id'].isin(chrom_bin_ids) &
            pixels['bin2_id'].isin(chrom_bin_ids)
        ].copy()

        if len(chrom_pixels) == 0:
            continue

        # Create mapping from global bin_id to genomic bin index (position / resolution)
        # scHiCluster expects indices based on genomic position, not sequential numbering
        bin_id_to_index = {bin_id: bins.loc[bin_id, 'start'] // resolution for bin_id in chrom_bin_ids}

        # Map bin IDs to genomic indices
        chrom_pixels['bin1_idx'] = chrom_pixels['bin1_id'].map(bin_id_to_index)
        chrom_pixels['bin2_idx'] = chrom_pixels['bin2_id'].map(bin_id_to_index)

        # Extract just the cell name without /cells/ prefix for filename
        cell_name_short = cell_name.split('/')[-1]

        # Write to text file in scHiCluster format: bin1_idx bin2_idx count
        output_file = os.path.join(output_dir, f"{cell_name_short}_{chrom}.txt")
        with open(output_file, 'w') as f:
            for _, row in chrom_pixels.iterrows():
                f.write(f"{int(row['bin1_idx'])}\t{int(row['bin2_idx'])}\t{int(row['count'])}\n")

        cells_exported += 1

    print(f"✅ Exported {cells_exported} cells for {chrom}")
    return cells_exported


def generate_chrom_sizes(scool_path, output_file):
    """
    Generate chromosome sizes file from scool.
    Only includes standard chromosomes (chr1-chr22/chr19, chrX, chrY).

    Args:
        scool_path: Path to scool file
        output_file: Path to output chromosome sizes file
    """
    import re

    # Open first cell to get chromosome info
    scool = list_scool_cells(scool_path)
    first_cell = scool[0]
    # cell path already includes /cells/ prefix
    clr = cooler.Cooler(f"{scool_path}::{first_cell}")

    # Pattern for standard chromosomes: chr1-chr99, chrX, chrY
    standard_chrom_pattern = re.compile(r'^chr(\d+|X|Y)$')

    with open(output_file, 'w') as f:
        for chrom, size in clr.chromsizes.items():
            # Only write standard chromosomes
            if standard_chrom_pattern.match(chrom):
                f.write(f"{chrom}\t{size}\n")

    print(f"✅ Created chromosome sizes file: {output_file} (standard chromosomes only)")


def main():
    parser = argparse.ArgumentParser(
        description="Export cells from scool to text format for scHiCluster"
    )
    parser.add_argument("--scool-file", required=True, help="Input scool file")
    parser.add_argument("--chrom", required=True, help="Chromosome to export")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--chrom-sizes-file",
        help="Generate chromosome sizes file (optional)"
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Maximum number of cells to export (for testing)"
    )

    args = parser.parse_args()

    try:
        # Generate chrom sizes if requested
        if args.chrom_sizes_file:
            generate_chrom_sizes(args.scool_file, args.chrom_sizes_file)

        # Export cells
        cells_exported = export_cells_for_chromosome(
            args.scool_file,
            args.chrom,
            args.output_dir,
            args.max_cells
        )

        if cells_exported == 0:
            print(f"⚠️  Warning: No cells exported for {args.chrom}")
            return 0  # Not an error - chromosome may have no data

        return 0

    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
