#!/usr/bin/env python
"""
Coarsen a scool file to lower resolution by aggregating bins.

This uses cooler's built-in coarsen functionality to change resolution
by a multiplier k (e.g., 10kb -> 60kb with k=6).

Usage:
    python coarsen_scool.py --input-scool <file.scool> --output-scool <coarsened.scool> --factor <k>
"""

import argparse
import cooler
import h5py
import sys
import os
import tempfile
from pathlib import Path


def coarsen_scool(input_scool, output_scool, factor):
    """
    Coarsen all cells in a scool file to lower resolution.

    Args:
        input_scool: Input scool file path
        output_scool: Output scool file path
        factor: Coarsening factor (k) - new_resolution = old_resolution * k

    Returns:
        Number of cells processed
    """
    # Get list of all cells
    try:
        cell_list = cooler.fileops.list_scool_cells(input_scool)
    except (OSError, AttributeError):
        # Fallback: read cells directly from HDF5
        with h5py.File(input_scool, 'r') as f:
            if 'cells' in f:
                cell_list = [f'/cells/{name}' for name in f['cells'].keys()]
            else:
                raise ValueError("No cells group found in scool file")

    n_cells = len(cell_list)

    print(f"Found {n_cells} cells in scool file")

    # Get original resolution
    first_cell = cell_list[0]
    clr = cooler.Cooler(f"{input_scool}::{first_cell}")
    original_resolution = clr.binsize
    new_resolution = original_resolution * factor

    print(f"Original resolution: {original_resolution}")
    print(f"New resolution: {new_resolution} (factor={factor})")
    print("")

    # Create temporary directory for coarsened cooler files
    temp_dir = tempfile.mkdtemp(prefix="coarsen_scool_")
    temp_coolers = []

    try:
        # Coarsen each cell
        for i, cell_name in enumerate(cell_list, 1):
            if i % 10 == 0 or i == 1 or i == n_cells:
                print(f"  Processing cell {i}/{n_cells}: {cell_name}")

            # Load cell cooler
            cell_uri = f"{input_scool}::{cell_name}"
            clr = cooler.Cooler(cell_uri)

            # Create temporary cooler file for this cell
            # Extract just the basename without the /cells/ prefix
            cell_basename = os.path.basename(cell_name)
            temp_cool = os.path.join(temp_dir, f"{cell_basename}.cool")

            # Coarsen the cooler
            cooler.coarsen_cooler(
                cell_uri,
                temp_cool,
                factor,
                chunksize=int(1e6)
            )

            temp_coolers.append((cell_name, temp_cool))

        print("")
        print(f"✅ Coarsened {n_cells} cells")
        print("")

        # Create new scool file
        print("Creating coarsened scool file...")

        # Create scool file manually using h5py
        with h5py.File(output_scool, 'w') as scool_file:
            # Create cells group
            cells_group = scool_file.create_group('cells')

            # Copy each cooler into the cells group
            for cell_name, temp_cool in temp_coolers:
                cell_basename = os.path.basename(cell_name)
                with h5py.File(temp_cool, 'r') as cool_file:
                    cool_file.copy(cool_file['/'], cells_group, name=cell_basename)

            # Copy bins and chroms from first cooler to root level
            first_cool = temp_coolers[0][1]
            with h5py.File(first_cool, 'r') as cool_file:
                cool_file.copy('bins', scool_file)
                cool_file.copy('chroms', scool_file)

        print(f"✅ Created: {output_scool}")

        # Verify output
        with h5py.File(output_scool, 'r') as f:
            if 'cells' in f:
                n_cells_out = len(list(f['cells'].keys()))
                print("")
                print("Verification:")
                print(f"  Cells in output: {n_cells_out}")

                # Check resolution from first cell
                first_cell_name = list(f['cells'].keys())[0]
                output_clr = cooler.Cooler(f"{output_scool}::cells/{first_cell_name}")
                output_resolution = output_clr.binsize
                print(f"  Output resolution: {output_resolution}")

        if output_resolution != new_resolution:
            print(f"  ⚠️  Warning: Expected {new_resolution}, got {output_resolution}")

        return n_cells

    finally:
        # Clean up temporary files
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Coarsen scool file to lower resolution"
    )
    parser.add_argument("--input-scool", required=True, help="Input scool file")
    parser.add_argument("--output-scool", required=True, help="Output coarsened scool file")
    parser.add_argument("--factor", type=int, required=True,
                       help="Coarsening factor k (new_res = old_res * k)")

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.input_scool):
        print(f"❌ Error: Input file not found: {args.input_scool}", file=sys.stderr)
        return 1

    if args.factor < 1:
        print(f"❌ Error: Factor must be >= 1, got {args.factor}", file=sys.stderr)
        return 1

    if args.factor == 1:
        print("⚠️  Warning: Factor=1 means no coarsening (copy would be created)")

    # Check if output already exists
    if os.path.exists(args.output_scool):
        print(f"⚠️  Warning: Output file already exists: {args.output_scool}")
        print("   It will be overwritten.")
        print("")

    try:
        print("=" * 80)
        print("Scool Coarsening")
        print("=" * 80)
        print(f"Input:  {args.input_scool}")
        print(f"Output: {args.output_scool}")
        print(f"Factor: {args.factor}")
        print("=" * 80)
        print("")

        n_cells = coarsen_scool(
            args.input_scool,
            args.output_scool,
            args.factor
        )

        print("")
        print("=" * 80)
        print("Coarsening Complete!")
        print("=" * 80)
        print(f"Processed {n_cells} cells")
        print(f"Output: {args.output_scool}")
        print("=" * 80)

        return 0

    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
