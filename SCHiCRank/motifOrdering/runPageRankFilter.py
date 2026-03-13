import argparse
import os
from runSCHiCRank import run_pagerank_filter


def main():
    parser = argparse.ArgumentParser(
        description="Run PageRank-based iterative cell filtering on pairwise similarity data."
    )
    parser.add_argument(
        "--input-dir", type=str, required=True,
        help="Directory containing pairwise similarity CSV files (e.g., 'pairwiseSimilarities/K4-long/')."
    )
    parser.add_argument(
        "--label", type=str, default=None,
        help="Label for output file naming (default: same as input directory name)."
    )
    parser.add_argument(
        "--k", type=int, default=5,
        help="Number of top neighbors to consider for each cell (default: 5)."
    )
    parser.add_argument(
        "--min-active-cells", type=int, default=10,
        help="Minimum number of active cells to keep in the analysis (default: 10)."
    )
    parser.add_argument(
        "--plots", dest="plots", action="store_true",
        help="Generate and display plots (default: enabled)."
    )
    parser.add_argument(
        "--no-plots", dest="plots", action="store_false",
        help="Disable plot generation."
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all output files (default: current directory)."
    )
    parser.set_defaults(plots=True)

    args = parser.parse_args()

    # Use input directory name as label if not specified
    label = args.label if args.label else args.input_dir.replace("/", "_").replace("\\", "_").strip("_")

    # Update global configuration in runSCHiCRank module
    import runSCHiCRank
    runSCHiCRank.K = args.k
    runSCHiCRank.MIN_ACTIVE_CELLS = args.min_active_cells

    # Set the metadata file path to work_dir
    if args.work_dir:
        runSCHiCRank.metaFn = os.path.join(args.work_dir, "cellAndPhaseInfo.pkl")

    print(f"Running PageRank filter on: {args.input_dir}")
    print(f"Configuration: K={args.k}, MIN_ACTIVE_CELLS={args.min_active_cells}, plots={args.plots}")

    run_pagerank_filter(args.input_dir, label=label, plots=args.plots, work_dir=args.work_dir)


if __name__ == "__main__":
    main()
