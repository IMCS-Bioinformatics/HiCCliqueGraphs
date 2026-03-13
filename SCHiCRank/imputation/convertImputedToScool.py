#!/usr/bin/env python
"""
Convert imputed HDF5 files from scHiCluster to scool format.

Usage:
    python convertImputedToScool.py --imputed-dir <dir> --output-scool <path> --chrom-sizes <file> --resolution <int>
"""

import argparse
import cooler
import h5py
import numpy as np
import pandas as pd
import os
import sys
import glob
from scipy.sparse import csr_matrix


def create_bins_dataframe(chrom_sizes_file, resolution):
    """
    Create bins DataFrame for all chromosomes.

    Args:
        chrom_sizes_file: Path to chromosome sizes file
        resolution: Bin resolution in base pairs

    Returns:
        DataFrame with columns: chrom, start, end
    """
    # Read chromosome sizes
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)

    # Create bins for all chromosomes
    all_bins = []
    for chrom in sorted(chrom_sizes.keys()):
        chrom_size = chrom_sizes[chrom]
        n_bins = int(np.ceil(chrom_size / resolution))

        chrom_bins = pd.DataFrame({
            'chrom': [chrom] * n_bins,
            'start': np.arange(n_bins) * resolution,
            'end': np.minimum(np.arange(1, n_bins + 1) * resolution, chrom_size)
        })
        all_bins.append(chrom_bins)

    bins_df = pd.concat(all_bins, ignore_index=True)
    return bins_df, chrom_sizes


def read_imputed_hdf5(hdf5_file):
    """
    Read imputed data from HDF5 file.

    Args:
        hdf5_file: Path to HDF5 file

    Returns:
        Sparse COO matrix
    """
    with h5py.File(hdf5_file, "r") as f:
        group = f["Matrix"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]

    matrix_sparse = csr_matrix((data, indices, indptr))
    return matrix_sparse.tocoo()


def process_imputed_cells(imputed_dir, bins_df, chromosomes):
    """
    Process all imputed HDF5 files and organize by cell.

    Args:
        imputed_dir: Directory containing imputed HDF5 files
        bins_df: DataFrame of bins
        chromosomes: List of chromosomes to process

    Returns:
        Dictionary mapping cell_name -> list of pixel DataFrames
    """
    cell_coolers = {}

    for chrom in chromosomes:
        chrom_dir = os.path.join(imputed_dir, chrom)
        if not os.path.exists(chrom_dir):
            print(f"⚠️  Warning: Directory not found: {chrom_dir}")
            continue

        hdf5_files = glob.glob(os.path.join(chrom_dir, "*.hdf5"))
        print(f"Found {len(hdf5_files)} imputed cells for {chrom}")

        # Get chromosome bin offset
        chrom_start_bin = bins_df[bins_df['chrom'] == chrom].index[0]

        for hdf5_file in hdf5_files:
            # Extract cell name from filename
            basename = os.path.basename(hdf5_file)
            # Typical format: cellname_chrX_params.hdf5
            parts = basename.replace('.hdf5', '').split('_')
            # Find where chromosome starts
            chrom_idx = next((i for i, p in enumerate(parts) if p.startswith('chr')), None)
            if chrom_idx is None:
                print(f"⚠️  Warning: Could not parse cell name from {basename}")
                continue

            cell_name = '_'.join(parts[:chrom_idx])

            if cell_name not in cell_coolers:
                cell_coolers[cell_name] = []

            # Read imputed data
            matrix_coo = read_imputed_hdf5(hdf5_file)

            # Create pixels with global bin IDs
            bin1 = matrix_coo.row + chrom_start_bin
            bin2 = matrix_coo.col + chrom_start_bin

            # Ensure bin1_id <= bin2_id (swap if needed)
            bin1_final = np.minimum(bin1, bin2)
            bin2_final = np.maximum(bin1, bin2)

            pixels = pd.DataFrame({
                'bin1_id': bin1_final,
                'bin2_id': bin2_final,
                'count': matrix_coo.data
            })

            cell_coolers[cell_name].append(pixels)

    return cell_coolers


def create_scool_from_cells(cell_coolers, bins_df, output_scool):
    """
    Create scool file from cell coolers.

    Args:
        cell_coolers: Dictionary mapping cell_name -> list of pixel DataFrames
        bins_df: Bins DataFrame
        output_scool: Path to output scool file
    """
    print(f"Combining {len(cell_coolers)} cells into scool...")

    temp_coolers = []

    for i, (cell_name, pixel_list) in enumerate(cell_coolers.items()):
        if i % 10 == 0:
            print(f"  Processing cell {i+1}/{len(cell_coolers)}: {cell_name}")

        # Combine all chromosomes
        all_pixels = pd.concat(pixel_list, ignore_index=True)

        # Group by bin pairs and sum counts (in case of duplicates from matrix symmetry)
        all_pixels = all_pixels.groupby(['bin1_id', 'bin2_id'], as_index=False)['count'].sum()
        all_pixels = all_pixels.sort_values(['bin1_id', 'bin2_id'])

        # Create temporary cooler
        temp_cool = f"{output_scool}.tmp.{i}.cool"
        cooler.create_cooler(
            temp_cool,
            bins=bins_df,
            pixels=all_pixels,
            dtypes={'count': 'float'},
            ordered=True
        )
        temp_coolers.append((cell_name, temp_cool))

    # Create scool file from individual coolers
    print("Creating final scool file...")
    # For scool, we need to copy each cooler into /cells/ group
    import h5py
    from datetime import datetime

    # Get bin info from first cooler
    first_cool = temp_coolers[0][1]
    with h5py.File(first_cool, 'r') as f:
        nbins = len(bins_df)
        nchroms = len(bins_df['chrom'].unique())
        bin_size = bins_df['end'].iloc[0] - bins_df['start'].iloc[0]

    # Create the scool file with proper metadata
    with h5py.File(output_scool, 'w') as scool_file:
        # Add scool format attributes
        scool_file.attrs['format'] = 'HDF5::SCOOL'
        scool_file.attrs['format-version'] = '1'
        scool_file.attrs['format-url'] = 'https://github.com/mirnylab/cooler'
        scool_file.attrs['bin-type'] = 'fixed'
        scool_file.attrs['bin-size'] = bin_size
        scool_file.attrs['ncells'] = len(temp_coolers)
        scool_file.attrs['nbins'] = nbins
        scool_file.attrs['nchroms'] = nchroms
        scool_file.attrs['creation-date'] = datetime.now().isoformat()
        scool_file.attrs['generated-by'] = 'convertImputedToScool.py'
        scool_file.attrs['genome-assembly'] = 'unknown'
        scool_file.attrs['metadata'] = '{}'

        # Copy bins and chroms from first cooler to root level
        # This is required for valid scool format
        with h5py.File(first_cool, 'r') as cool_file:
            cool_file.copy('bins', scool_file)
            cool_file.copy('chroms', scool_file)

        # Create cells group and copy all cells
        cells_group = scool_file.create_group('cells')

        for cell_name, temp_cool_path in temp_coolers:
            # Copy the entire cooler structure into /cells/{cell_name}
            with h5py.File(temp_cool_path, 'r') as cool_file:
                cool_file.copy(cool_file['/'], cells_group, name=cell_name)

    # Clean up temporary files
    print("Cleaning up temporary files...")
    for cell, temp_cool in temp_coolers:
        if os.path.exists(temp_cool):
            os.remove(temp_cool)

    print(f"✅ Created imputed scool: {output_scool}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert imputed HDF5 files to scool format"
    )
    parser.add_argument("--imputed-dir", required=True, help="Directory with imputed HDF5 files")
    parser.add_argument("--output-scool", required=True, help="Output scool file path")
    parser.add_argument("--chrom-sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--resolution", type=int, required=True, help="Bin resolution in bp")
    parser.add_argument("--chromosomes", nargs='+', required=True, help="List of chromosomes")

    args = parser.parse_args()

    try:
        # Create bins DataFrame
        print("Creating bins...")
        bins_df, chrom_sizes = create_bins_dataframe(args.chrom_sizes, args.resolution)
        print(f"Created {len(bins_df)} bins across {len(chrom_sizes)} chromosomes")

        # Process imputed cells
        print("Processing imputed cells...")
        cell_coolers = process_imputed_cells(args.imputed_dir, bins_df, args.chromosomes)

        if not cell_coolers:
            print("❌ No cells found to process")
            return 1

        # Create scool
        create_scool_from_cells(cell_coolers, bins_df, args.output_scool)

        return 0

    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
