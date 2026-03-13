import argparse
from processOriginalCoolDataset import process_cells

ALL_CHROMOSOMES = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chrX",
]


def main():
    parser = argparse.ArgumentParser(
        description="Process a multicool scHi-C dataset into per-chromosome interaction files."
    )
    parser.add_argument(
        "--cell-count", type=int, default=None,
        help="Number of cells to process (default: all cells)."
    )
    parser.add_argument(
        "-k", type=int, default=6,
        help="Resolution multiplier applied to the base resolution (default: 6)."
    )
    parser.add_argument(
        "--chromosomes", nargs="+", default=ALL_CHROMOSOMES,
        help="Chromosomes to process. Use 'ALL' for all chromosomes (default: ALL)."
    )
    parser.add_argument(
        "--fn", type=str, default="sourceData/nagano_10kb_cell_types.scool",
        help="Path to the input .scool file (default: sourceData/nagano_10kb_cell_types.scool)."
    )
    parser.add_argument(
        "--fn-resolution", type=int, default=10000,
        help="Base resolution of the .scool file in bp (default: 10000)."
    )
    parser.add_argument(
        "--postfix", type=str, default="base60k",
        help="Postfix used in output file names (default: base60k)."
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all output files (default: current directory)."
    )

    args = parser.parse_args()

    # Handle "ALL" keyword for chromosomes
    chromosomes = ALL_CHROMOSOMES if args.chromosomes == ["ALL"] else args.chromosomes

    process_cells(
        cellCount=args.cell_count,
        k=args.k,
        chromosomes=chromosomes,
        fn=args.fn,
        fnResolution=args.fn_resolution,
        postfix=args.postfix,
        work_dir=args.work_dir,
    )


if __name__ == "__main__":
    main()
