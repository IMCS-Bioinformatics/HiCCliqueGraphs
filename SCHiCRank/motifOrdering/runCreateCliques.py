import argparse
import os
from createCliqueDatafiles import createCliquePickles

ALL_CHROMOSOMES = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chrX",
]


def main():
    parser = argparse.ArgumentParser(
        description="Generate clique data (K3-K8) for processed scHi-C chromosome files."
    )
    parser.add_argument(
        "--chromosomes", nargs="+", default=ALL_CHROMOSOMES,
        help="Chromosomes to process. Use 'ALL' for all chromosomes (default: ALL)."
    )
    parser.add_argument(
        "--postfix", type=str, default="base60k",
        help="Postfix used in input file names (default: base60k)."
    )
    parser.add_argument(
        "--resolution", type=int, default=60000,
        help="Resolution in bp, used in file naming (default: 60000)."
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all input/output files (default: current directory)."
    )

    args = parser.parse_args()

    # Handle "ALL" keyword for chromosomes
    chromosomes = ALL_CHROMOSOMES if args.chromosomes == ["ALL"] else args.chromosomes

    for ch in chromosomes:
        baseFn = os.path.join(args.work_dir, f"{args.postfix}-{ch}-{args.resolution}.pkl")
        resFn = os.path.join(args.work_dir, f"{args.postfix}-{ch}-{args.resolution}-cliques.pkl")
        print(f"Processing {ch}...")
        createCliquePickles(baseFn, resFn)

    print("All chromosomes processed.")


if __name__ == "__main__":
    main()
