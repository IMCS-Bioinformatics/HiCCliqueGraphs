import argparse
import pandas as pd
import numpy as np
from itertools import permutations
from scipy.stats import kendalltau
from datetime import datetime
import pickle
import os


def load_cell_phases(meta_file="cellAndPhaseInfo.pkl"):
    """Load cell phase metadata."""
    with open(meta_file, "rb") as f:
        meta = pickle.load(f)
    return meta["cell_phase"]


def collect_motif_counts_per_chromosome(chromosomes, postfix, resolution, work_dir=""):
    """
    Collect motif counts from clique pickle files for all chromosomes.

    Returns:
        dict: {chromosome: {motif_type: {cell_id: count}}}
    """
    all_data = {}

    for chrom in chromosomes:
        clique_file = os.path.join(work_dir, f"{postfix}-{chrom}-{resolution}-cliques.pkl")

        if not os.path.exists(clique_file):
            print(f"Warning: {clique_file} not found, skipping...")
            continue

        print(f"Reading {clique_file}...")

        with open(clique_file, "rb") as f:
            data = pickle.load(f)

        chrom_data = {}

        # Process each K type (K3, K4, ..., K8)
        for k_type in data["cell_cliques"].keys():
            # Process both alllengths and long
            for length_type in ["alllengths", "long"]:
                motif_name = f"{k_type}-{length_type}"
                cell_counts = {}

                # Get cliques for each cell
                for cell_id, cliques in data["cell_cliques"][k_type].items():
                    if length_type == "long":
                        # Count only long cliques (span >= 2Mb)
                        count = len([c for c in cliques if c[-1] - c[0] >= 2000000])
                    else:
                        # Count all cliques
                        count = len(cliques)

                    cell_counts[cell_id] = count

                chrom_data[motif_name] = cell_counts

        all_data[chrom] = chrom_data

    return all_data


def rank_cells_by_motif_count(motif_counts_per_chrom, motif_type, min_count=10):
    """
    Rank cells by motif count, aggregating across chromosomes using median rank.
    """
    cell_ranks = {}  # {cell_id: [ranks across chromosomes]}

    # For each chromosome, rank cells by motif count
    for chrom, motif_data in motif_counts_per_chrom.items():
        if motif_type not in motif_data:
            continue

        counts = motif_data[motif_type]

        # Filter cells with sufficient motif count
        eligible_cells = {cell_id: count for cell_id, count in counts.items()
                         if count >= min_count}

        if not eligible_cells:
            continue

        # Sort by count (descending) and assign ranks
        sorted_cells = sorted(eligible_cells.items(), key=lambda x: x[1], reverse=True)

        for rank, (cell_id, count) in enumerate(sorted_cells):
            if cell_id not in cell_ranks:
                cell_ranks[cell_id] = []
            cell_ranks[cell_id].append(rank)

    # Calculate median rank for each cell
    cell_median_ranks = []
    for cell_id, ranks in cell_ranks.items():
        if len(ranks) > 0:  # Only include cells that appeared in at least one chromosome
            median_rank = np.median(ranks)
            cell_median_ranks.append((cell_id, median_rank))

    # Sort by median rank
    cell_median_ranks.sort(key=lambda x: x[1])

    return cell_median_ranks


def calculate_tau_for_phase_ordering(ranked_cells, cell_phases, phase_order):
    """Calculate Kendall's tau for a specific phase ordering."""
    phase_to_rank = {phase: i for i, phase in enumerate(phase_order)}

    # Create ideal ranks based on phase ordering
    ideal_ranks = []
    actual_ranks = []

    for i, (cell_id, _) in enumerate(ranked_cells):
        phase = cell_phases[cell_id] if cell_id < len(cell_phases) else "Unknown"
        ideal_rank = phase_to_rank.get(phase, len(phase_order))

        ideal_ranks.append(ideal_rank)
        actual_ranks.append(i)

    tau, p_value = kendalltau(actual_ranks, ideal_ranks)
    return tau, p_value


def calculate_best_tau(ranked_cells, cell_phases):
    """Try all phase orderings and return the best tau."""
    unique_phases = set(cell_phases[cell_id] for cell_id, _ in ranked_cells
                       if cell_id < len(cell_phases))

    best_tau = -2
    best_ordering = None
    best_p_value = None

    for phase_order in permutations(unique_phases):
        tau, p_value = calculate_tau_for_phase_ordering(ranked_cells, cell_phases, phase_order)
        if tau > best_tau:
            best_tau = tau
            best_ordering = phase_order
            best_p_value = p_value

    return best_tau, best_ordering, best_p_value


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Kendall's tau for motif-count-based cell orderings."
    )
    parser.add_argument(
        "--motif-type", type=str, required=True,
        help="Motif type to analyze (e.g., 'K3-alllengths', 'K4-long')."
    )
    parser.add_argument(
        "--postfix", type=str, default="base60k",
        help="Postfix used in clique file names (default: base60k)."
    )
    parser.add_argument(
        "--resolution", type=int, default=60000,
        help="Resolution in bp (default: 60000)."
    )
    parser.add_argument(
        "--min-count", type=int, default=10,
        help="Minimum motif count per chromosome to include cell (default: 10)."
    )
    parser.add_argument(
        "--log-file", type=str, default="tau_scores.txt",
        help="Log file to append results (default: tau_scores.txt)."
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all input/output files (default: current directory)."
    )

    args = parser.parse_args()

    ALL_CHROMOSOMES = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chrX",
    ]

    meta_file = os.path.join(args.work_dir, "cellAndPhaseInfo.pkl") if args.work_dir else "cellAndPhaseInfo.pkl"
    log_file = os.path.join(args.work_dir, args.log_file) if args.work_dir else args.log_file

    print("="*70)
    print(f"MOTIF-COUNT-BASED ORDERING: {args.motif_type}")
    print("="*70)

    # Load cell phases
    print("Loading cell phase metadata...")
    cell_phases = load_cell_phases(meta_file)

    # Collect motif counts
    print("Collecting motif counts from all chromosomes...")
    motif_counts = collect_motif_counts_per_chromosome(
        ALL_CHROMOSOMES, args.postfix, args.resolution, work_dir=args.work_dir
    )

    # Rank cells by motif count
    print(f"Ranking cells by {args.motif_type} count (min_count={args.min_count})...")
    ranked_cells = rank_cells_by_motif_count(motif_counts, args.motif_type, args.min_count)

    print(f"Ranked {len(ranked_cells)} cells")

    if len(ranked_cells) == 0:
        print("ERROR: No cells to rank!")
        return

    if len(ranked_cells) < 3:
        print(f"WARNING: Only {len(ranked_cells)} ranked cells - too few for Kendall's tau (need >= 3). Skipping.")
        return

    # Calculate best tau
    print("Calculating Kendall's tau for all phase orderings...")
    best_tau, best_ordering, best_p_value = calculate_best_tau(ranked_cells, cell_phases)

    if best_ordering is None or best_p_value is None:
        print("WARNING: Could not compute Kendall's tau (not enough unique phases or cells). Skipping.")
        return

    print("\n" + "="*70)
    print("RESULT:")
    print(f"Best Kendall's tau: {best_tau:.4f} (p-value: {best_p_value:.4e})")
    print(f"Best phase ordering: {' -> '.join(best_ordering)}")
    print("="*70)

    # Log to file
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"""
{'='*70}
Timestamp: {timestamp}
Method: Motif-Count-Based Ordering
Configuration:
  - Motif: {args.motif_type}
  - Postfix: {args.postfix}
  - Resolution: {args.resolution}
  - Min count threshold: {args.min_count}

Results:
  - Number of ranked cells: {len(ranked_cells)}
  - Best Kendall's tau: {best_tau:.4f} (p-value: {best_p_value:.4e})
  - Best phase ordering: {' -> '.join(best_ordering)}
{'='*70}

"""

    with open(log_file, "a") as f:
        f.write(log_entry)

    print(f"\nResults appended to {log_file}")


if __name__ == "__main__":
    main()
