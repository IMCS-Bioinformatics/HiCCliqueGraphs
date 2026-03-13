import argparse
import os
import pandas as pd
from itertools import permutations
from scipy.stats import kendalltau
from datetime import datetime


def calculate_tau_for_phase_ordering(df, phase_order):
    """
    Calculate Kendall's tau for a specific phase ordering.

    Args:
        df: DataFrame with 'Phase' column (already sorted by actual ranking)
        phase_order: Tuple of phases in desired order

    Returns:
        Kendall's tau correlation coefficient
    """
    # Create ideal ranking where cells are ordered by phase_order
    phase_to_rank = {phase: i for i, phase in enumerate(phase_order)}

    # Assign ideal ranks to each cell based on its phase
    ideal_ranks = [phase_to_rank.get(phase, len(phase_order)) for phase in df['Phase']]

    # Actual ranks are just the row indices (0, 1, 2, ...)
    actual_ranks = list(range(len(df)))

    # Calculate Kendall's tau
    tau, p_value = kendalltau(actual_ranks, ideal_ranks)

    return tau, p_value


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Kendall's tau to evaluate phase ordering quality in cell rankings."
    )
    parser.add_argument(
        "--input-csv", type=str, required=True,
        help="Path to the CSV file from step 4 (e.g., final_active_cells_*.csv)."
    )
    parser.add_argument(
        "--motif-config", type=str, default="unknown",
        help="Motif configuration used in step 4 (e.g., 'K3-alllengths')."
    )
    parser.add_argument(
        "--k", type=int, default=5,
        help="K value used in step 4 (number of neighbors)."
    )
    parser.add_argument(
        "--min-active-cells", type=int, default=10,
        help="Min active cells value used in step 4."
    )
    parser.add_argument(
        "--log-file", type=str, default="tau_scores.txt",
        help="Path to log file for results (default: tau_scores.txt)."
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all input/output files (default: current directory)."
    )

    args = parser.parse_args()

    input_csv = os.path.join(args.work_dir, args.input_csv) if args.work_dir else args.input_csv
    log_file = os.path.join(args.work_dir, args.log_file) if args.work_dir else args.log_file

    # Read the CSV file
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)

    # Get unique phases
    unique_phases = df['Phase'].unique()
    print(f"Found {len(unique_phases)} unique phases: {list(unique_phases)}")

    # Try all permutations of phase orderings
    print(f"Testing all {len(list(permutations(unique_phases)))} possible phase orderings...")

    best_tau = -2  # Start with impossible value
    best_ordering = None
    best_p_value = None
    all_results = []

    for phase_order in permutations(unique_phases):
        tau, p_value = calculate_tau_for_phase_ordering(df, phase_order)
        all_results.append((tau, p_value, phase_order))

        if tau > best_tau:
            best_tau = tau
            best_ordering = phase_order
            best_p_value = p_value

    # Sort results by tau (descending)
    all_results.sort(key=lambda x: x[0], reverse=True)

    # Print results
    print("\n" + "="*70)
    print("BEST RESULT:")
    print(f"Kendall's tau: {best_tau:.4f} (p-value: {best_p_value:.4e})")
    print(f"Best phase ordering: {' -> '.join(best_ordering)}")
    print("="*70)

    print("\nTop 5 phase orderings:")
    for i, (tau, p_val, ordering) in enumerate(all_results[:5], 1):
        print(f"{i}. tau={tau:.4f} (p={p_val:.4e}): {' -> '.join(ordering)}")

    # Log to file
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"""
{'='*70}
Timestamp: {timestamp}
Input CSV: {args.input_csv}
Configuration:
  - Motif: {args.motif_config}
  - K (neighbors): {args.k}
  - Min active cells: {args.min_active_cells}

Results:
  - Number of cells: {len(df)}
  - Number of phases: {len(unique_phases)}
  - Phases found: {', '.join(unique_phases)}

Best Kendall's tau: {best_tau:.4f} (p-value: {best_p_value:.4e})
Best phase ordering: {' -> '.join(best_ordering)}

Top 5 orderings:
"""
    for i, (tau, p_val, ordering) in enumerate(all_results[:5], 1):
        log_entry += f"  {i}. tau={tau:.4f} (p={p_val:.4e}): {' -> '.join(ordering)}\n"

    log_entry += "="*70 + "\n"

    # Append to log file
    with open(log_file, "a") as f:
        f.write(log_entry)

    print(f"\nResults appended to {log_file}")


if __name__ == "__main__":
    main()
