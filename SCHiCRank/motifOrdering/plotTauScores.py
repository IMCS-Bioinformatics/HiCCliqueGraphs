import re
import os
import matplotlib.pyplot as plt
import argparse


def parse_tau_scores(log_file):
    """
    Parse tau_scores.txt and extract motif configurations and their best tau scores.

    Returns:
        tuple: (dict of {motif_config: best_tau_score}, baseline_tau or None)
    """
    results = {}
    baseline_tau = None

    with open(log_file, "r") as f:
        content = f.read()

    # Split by the separator lines
    entries = content.split("=" * 70)

    for entry in entries:
        if not entry.strip():
            continue

        # Check if this is the baseline entry
        if "RANDOM BASELINE CALCULATION" in entry:
            baseline_match = re.search(r"Mean tau \(random baseline\):\s*([-+]?\d*\.\d+)", entry)
            if baseline_match:
                baseline_tau = float(baseline_match.group(1))
                continue

        # Check if this is motif-count-based or PageRank-based
        is_motif_count = "Motif-Count-Based Ordering" in entry

        # Extract motif configuration
        motif_match = re.search(r"- Motif:\s*(\S+)", entry)
        # Extract best tau score
        tau_match = re.search(r"Best Kendall's tau:\s*([-+]?\d*\.\d+)", entry)

        if motif_match and tau_match:
            motif = motif_match.group(1)
            tau = float(tau_match.group(1))

            # Add suffix to distinguish motif-count-based from PageRank-based
            if is_motif_count:
                motif = f"{motif}\n(count)"

            results[motif] = tau

    return results, baseline_tau


def plot_tau_scores(results, baseline_tau=None, output_file=None):
    """
    Create a bar chart of tau scores by motif configuration.
    """
    # Sort configurations for better visualization
    def sort_key(item):
        motif = item[0]
        base_motif = motif.replace("\n(count)", "")
        prefix = base_motif.split("-")[0]
        if prefix.startswith("K") and prefix[1:].isdigit():
            k_num = int(prefix[1:])
        elif prefix == "interaction":
            k_num = 100  # sort interactions after all K motifs
        else:
            k_num = 999
        length_type = 0 if "alllengths" in base_motif else 1
        method_type = 1 if "(count)" in motif else 0
        return (k_num, length_type, method_type)

    sorted_results = sorted(results.items(), key=sort_key)

    # Add baseline at the beginning if available
    if baseline_tau is not None:
        motifs = ["Random\nBaseline"] + [item[0] for item in sorted_results]
        taus = [baseline_tau] + [item[1] for item in sorted_results]
    else:
        motifs = [item[0] for item in sorted_results]
        taus = [item[1] for item in sorted_results]

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create bar chart
    bars = ax.bar(range(len(motifs)), taus, color='steelblue', alpha=0.8)

    # Color bars
    for i, motif in enumerate(motifs):
        if motif == "Random\nBaseline":
            bars[i].set_color('gray')
            bars[i].set_alpha(0.6)
        elif "(count)" in motif:
            if "alllengths" in motif:
                bars[i].set_color('lightblue')
                bars[i].set_edgecolor('steelblue')
                bars[i].set_linewidth(1.5)
            else:
                bars[i].set_color('lightsalmon')
                bars[i].set_edgecolor('coral')
                bars[i].set_linewidth(1.5)
        elif "alllengths" in motif:
            bars[i].set_color('steelblue')
        else:
            bars[i].set_color('coral')

    # Customize plot
    ax.set_xlabel("Configuration", fontsize=12, fontweight='bold')
    ax.set_ylabel("Best Kendall's Tau Score", fontsize=12, fontweight='bold')
    ax.set_title("Cell Phase Ordering Quality by Motif Type", fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(motifs)))
    ax.set_xticklabels(motifs, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Add legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor='steelblue', alpha=0.8, label='PageRank: All lengths'),
        Patch(facecolor='coral', alpha=0.8, label='PageRank: Long (>2Mb)'),
        Patch(facecolor='lightblue', edgecolor='steelblue', linewidth=1.5, label='Count: All lengths'),
        Patch(facecolor='lightsalmon', edgecolor='coral', linewidth=1.5, label='Count: Long (>2Mb)')
    ]
    if baseline_tau is not None:
        legend_elements.insert(0, Patch(facecolor='gray', alpha=0.6, label='Random baseline'))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # Add value labels on bars
    for i, (motif, tau) in enumerate(zip(motifs, taus)):
        ax.text(i, tau + 0.005, f'{tau:.3f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Visualize Kendall's tau scores from tau_scores.txt"
    )
    parser.add_argument(
        "--input", type=str, default="tau_scores.txt",
        help="Path to tau scores log file (default: tau_scores.txt)"
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Path to save figure (optional; if not provided, shows plot interactively)"
    )
    parser.add_argument(
        "--work-dir", type=str, default="",
        help="Directory for all input/output files (default: current directory)."
    )

    args = parser.parse_args()

    input_file = os.path.join(args.work_dir, args.input) if args.work_dir else args.input
    output_file = os.path.join(args.work_dir, args.output) if (args.work_dir and args.output) else args.output

    # Parse results
    print(f"Reading {input_file}...")
    results, baseline_tau = parse_tau_scores(input_file)

    if not results:
        print("No results found in the log file!")
        return

    if baseline_tau is not None:
        print(f"Random baseline tau: {baseline_tau:.4f}\n")

    print(f"Found {len(results)} configurations:")
    for motif, tau in sorted(results.items()):
        print(f"  {motif}: {tau:.4f}")

    # Create plot
    print("\nGenerating plot...")
    plot_tau_scores(results, baseline_tau, output_file)


if __name__ == "__main__":
    main()
