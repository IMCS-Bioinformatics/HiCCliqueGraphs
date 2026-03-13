import argparse
import os
import pickle
from MulticoolProcessor import MulticoolProcessor


def generate_cell_phase_info(
    scool_file="sourceData/nagano_10kb_cell_types.scool",
    cell_type_file="sourceData/nagano_assoziated_cell_types.txt",
    output_file="cellAndPhaseInfo.pkl"
):
    """
    Generate cellAndPhaseInfo.pkl from source data files.

    This script reads the .scool file to get cell names in order, then maps
    each cell to its cell cycle phase from the cell type annotation file.

    Args:
        scool_file: Path to the .scool file containing cell data
        cell_type_file: Path to tab-separated file mapping cell names to phases
        output_file: Path for the output pickle file
    """
    print("="*70)
    print("Generating Cell Phase Info")
    print("="*70)

    # Read cell names from .scool file
    print(f"\nReading cell names from {scool_file}...")
    M = MulticoolProcessor(scool_file)
    cell_names = M.getCellNames()
    print(f"Found {len(cell_names)} cells")

    # Read cell type mapping
    print(f"\nReading cell phase annotations from {cell_type_file}...")
    cell_type_mapping = {}

    with open(cell_type_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) == 2:
                    cell_name, cell_type = parts
                    # Remove the "/cells/" prefix if present
                    cell_name = cell_name.lstrip("/cells/")
                    cell_type_mapping[cell_name] = cell_type

    print(f"Found annotations for {len(cell_type_mapping)} cells")

    # Create ordered list of cell phases matching cell names order
    cell_phases = []
    missing_count = 0

    for cell_name in cell_names:
        if cell_name in cell_type_mapping:
            cell_phases.append(cell_type_mapping[cell_name])
        else:
            # If cell not found in annotations, mark as unknown
            cell_phases.append("unknown")
            missing_count += 1

    if missing_count > 0:
        print(f"\nWarning: {missing_count} cells have no phase annotation (marked as 'unknown')")

    # Count phases
    phase_counts = {}
    for phase in cell_phases:
        phase_counts[phase] = phase_counts.get(phase, 0) + 1

    print("\nCell cycle phase distribution:")
    for phase, count in sorted(phase_counts.items()):
        print(f"  {phase}: {count} cells ({100*count/len(cell_phases):.1f}%)")

    # Create the data structure
    data = {
        "cell_names": cell_names,
        "cell_phase": cell_phases
    }

    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Save to pickle file
    print(f"\nSaving to {output_file}...")
    with open(output_file, 'wb') as f:
        pickle.dump(data, f)

    print("\n" + "="*70)
    print("SUCCESS!")
    print(f"Created {output_file} with {len(cell_names)} cells")
    print("="*70)


def main():
    parser = argparse.ArgumentParser(
        description="Generate cellAndPhaseInfo.pkl from source data files."
    )
    parser.add_argument(
        "--scool-file",
        type=str,
        default="sourceData/nagano_10kb_cell_types.scool",
        help="Path to the .scool file (default: sourceData/nagano_10kb_cell_types.scool)"
    )
    parser.add_argument(
        "--cell-type-file",
        type=str,
        default="sourceData/nagano_assoziated_cell_types.txt",
        help="Path to the cell type annotation file (default: sourceData/nagano_assoziated_cell_types.txt)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="cellAndPhaseInfo.pkl",
        help="Output pickle file path (default: cellAndPhaseInfo.pkl)"
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for output files (default: current directory)."
    )

    args = parser.parse_args()

    output_file = args.output
    if args.work_dir:
        os.makedirs(args.work_dir, exist_ok=True)
        output_file = os.path.join(args.work_dir, args.output)

    generate_cell_phase_info(
        scool_file=args.scool_file,
        cell_type_file=args.cell_type_file,
        output_file=output_file
    )


if __name__ == "__main__":
    main()
