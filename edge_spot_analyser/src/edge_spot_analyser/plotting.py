"""
Figure plotting and statistical analysis module.

This module creates publication-ready figures from aggregated data:
- Small points: individual XY measurements
- Big points: replicate (date) means
- One-way ANOVA for statistical comparison between conditions
"""

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


@dataclass
class ColumnInfo:
    """Parsed information from a column name."""

    condition: str
    date: str
    well: str
    original: str


def parse_column_name(col: str) -> ColumnInfo | None:
    """
    Parse a column name into condition, date, and well.

    Expected format: {condition}_{date}_{well}
    e.g., "no_TRAK_230912_B02" -> condition="no_TRAK", date="230912", well="B02"

    Args:
        col: Column name string

    Returns:
        ColumnInfo or None if parsing fails
    """
    # Match pattern: anything_6digits[_suffix]_wellID (e.g., B02, G10)
    # Date can have optional suffix like 231120_1 for multiple experiments on same date
    match = re.match(r"^(.+)_(\d{6}(?:_\d+)?)_([A-H]\d{2})$", col)
    if not match:
        return None

    return ColumnInfo(
        condition=match.group(1),
        date=match.group(2),
        well=match.group(3),
        original=col,
    )


def restructure_data_for_plotting(
    df: pd.DataFrame,
    conditions_order: list[str] | None = None,
) -> dict[str, dict]:
    """
    Restructure figure data for plotting.

    Args:
        df: DataFrame with columns like {condition}_{date}_{well}
        conditions_order: Optional list of conditions in desired order

    Returns:
        Dict mapping condition -> {
            'individual_values': list of all XY values,
            'replicate_means': list of mean per date,
            'replicate_dates': list of date strings,
        }
    """
    # Parse all columns
    parsed_cols = {}
    for col in df.columns:
        if col == "XY":
            continue
        info = parse_column_name(col)
        if info:
            parsed_cols[col] = info

    # Group by condition, preserving order of first appearance in columns
    conditions = {}
    condition_order_seen = []
    for col, info in parsed_cols.items():
        if info.condition not in conditions:
            conditions[info.condition] = {"columns": [], "dates": set()}
            condition_order_seen.append(info.condition)
        conditions[info.condition]["columns"].append(col)
        conditions[info.condition]["dates"].add(info.date)

    # Order conditions
    if conditions_order:
        ordered_conditions = [c for c in conditions_order if c in conditions]
        # Add any conditions not in the order list
        ordered_conditions.extend([c for c in conditions if c not in ordered_conditions])
    else:
        # Preserve order from CSV columns (control is typically first)
        ordered_conditions = condition_order_seen

    # Build result
    result = {}
    for condition in ordered_conditions:
        cols = conditions[condition]["columns"]

        # Individual values: all non-NaN values from all columns
        individual_values = df[cols].values.flatten()
        individual_values = individual_values[~np.isnan(individual_values)]

        # Replicate means: mean per date
        dates_data = {}
        for col in cols:
            info = parsed_cols[col]
            if info.date not in dates_data:
                dates_data[info.date] = []
            # Get all non-NaN values for this well
            well_values = df[col].dropna().values
            dates_data[info.date].extend(well_values)

        replicate_means = []
        replicate_dates = []
        for date in sorted(dates_data.keys()):
            values = dates_data[date]
            if values:
                replicate_means.append(np.mean(values))
                replicate_dates.append(date)

        result[condition] = {
            "individual_values": individual_values.tolist(),
            "replicate_means": replicate_means,
            "replicate_dates": replicate_dates,
        }

    return result


def run_one_way_anova(data: dict[str, dict]) -> dict:
    """
    Run one-way ANOVA comparing conditions using replicate means.

    Args:
        data: Output from restructure_data_for_plotting

    Returns:
        Dict with ANOVA results:
            - f_statistic: F-statistic
            - p_value: p-value
            - conditions: list of conditions
            - n_per_condition: list of sample sizes (replicates)
            - means: list of grand means per condition
            - sems: list of SEM per condition
    """
    conditions = list(data.keys())
    groups = [data[c]["replicate_means"] for c in conditions]

    # Filter out empty groups
    valid_groups = [(c, g) for c, g in zip(conditions, groups, strict=False) if len(g) > 0]
    if len(valid_groups) < 2:
        logger.warning("Need at least 2 groups with data for ANOVA")
        return {
            "f_statistic": np.nan,
            "p_value": np.nan,
            "conditions": conditions,
            "n_per_condition": [len(g) for g in groups],
            "means": [np.mean(g) if g else np.nan for g in groups],
            "sems": [stats.sem(g) if len(g) > 1 else np.nan for g in groups],
        }

    valid_conditions, valid_data = zip(*valid_groups, strict=False)

    # Run one-way ANOVA
    f_stat, p_val = stats.f_oneway(*valid_data)

    return {
        "f_statistic": f_stat,
        "p_value": p_val,
        "conditions": conditions,
        "n_per_condition": [len(g) for g in groups],
        "means": [np.mean(g) if g else np.nan for g in groups],
        "sems": [stats.sem(g) if len(g) > 1 else np.nan for g in groups],
    }


def create_figure(
    data: dict[str, dict],
    title: str = "",
    ylabel: str = "Value",
    figsize: tuple[float, float] = (6, 5),
    small_point_size: float = 20,
    big_point_size: float = 100,
    small_point_alpha: float = 0.3,
    big_point_alpha: float = 0.8,
    jitter_width: float = 0.2,
    colors: list[str] | None = None,
):
    """
    Create a publication-ready figure with individual points and replicate means.

    Args:
        data: Output from restructure_data_for_plotting
        title: Figure title
        ylabel: Y-axis label
        figsize: Figure size in inches
        small_point_size: Size of individual XY points
        big_point_size: Size of replicate mean points
        small_point_alpha: Transparency of individual points
        big_point_alpha: Transparency of replicate mean points
        jitter_width: Width of horizontal jitter for points
        colors: Optional list of colors per condition

    Returns:
        matplotlib Figure object
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)

    conditions = list(data.keys())
    n_conditions = len(conditions)

    # Default colors
    if colors is None:
        colors = plt.cm.tab10.colors[:n_conditions]

    x_positions = np.arange(n_conditions)

    for i, (_condition, cond_data) in enumerate(data.items()):
        color = colors[i % len(colors)]

        # Plot individual values (small points) with jitter
        individual = cond_data["individual_values"]
        if individual:
            jitter = np.random.uniform(-jitter_width, jitter_width, len(individual))
            ax.scatter(
                x_positions[i] + jitter,
                individual,
                s=small_point_size,
                c=[color],
                alpha=small_point_alpha,
                edgecolors="none",
                zorder=1,
            )

        # Plot replicate means (big points) with less jitter
        replicate_means = cond_data["replicate_means"]
        if replicate_means:
            jitter = np.random.uniform(-jitter_width / 2, jitter_width / 2, len(replicate_means))
            ax.scatter(
                x_positions[i] + jitter,
                replicate_means,
                s=big_point_size,
                c=[color],
                alpha=big_point_alpha,
                edgecolors="black",
                linewidths=1,
                zorder=2,
            )

    # Configure axes
    ax.set_xticks(x_positions)
    ax.set_xticklabels(conditions, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Add horizontal line at y=1 for normalized data
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5, zorder=0)

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    return fig


def process_figure_file(
    csv_path: Path,
    output_dir: Path,
    conditions_order: list[str] | None = None,
    ylabel: str = "Value",
) -> dict:
    """
    Process a single figure CSV file and generate plot + statistics.

    Args:
        csv_path: Path to the figure CSV file (e.g., Fig1D_edge_spot_normalized.csv)
        output_dir: Output directory for plots and stats
        conditions_order: Optional list of conditions in desired order
        ylabel: Y-axis label for the plot

    Returns:
        Dict with ANOVA results
    """
    import matplotlib.pyplot as plt

    logger.info(f"Processing {csv_path.name}")

    # Load data
    df = pd.read_csv(csv_path, index_col=0)

    # Restructure for plotting
    data = restructure_data_for_plotting(df, conditions_order)

    if not data:
        logger.warning(f"No valid data found in {csv_path}")
        return {}

    # Run ANOVA
    anova_results = run_one_way_anova(data)

    # Extract figure name and metric from filename
    # e.g., "Fig1D_edge_spot_normalized.csv" -> fig_name="Fig1D", metric="edge_spot_normalized"
    stem = csv_path.stem
    parts = stem.split("_", 1)
    fig_name = parts[0]
    metric = parts[1] if len(parts) > 1 else stem

    # Create title
    title = f"{fig_name}: {metric.replace('_', ' ').title()}"
    if not np.isnan(anova_results["p_value"]):
        title += f"\nANOVA p={anova_results['p_value']:.4f}"

    # Create and save figure
    fig = create_figure(data, title=title, ylabel=ylabel)
    output_path = output_dir / f"{stem}_plot.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved plot to {output_path}")

    # Save PDF version
    output_path_pdf = output_dir / f"{stem}_plot.pdf"
    fig = create_figure(data, title=title, ylabel=ylabel)
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved PDF to {output_path_pdf}")

    # Save statistics
    stats_path = output_dir / f"{stem}_stats.csv"
    stats_df = pd.DataFrame(
        {
            "condition": anova_results["conditions"],
            "n_replicates": anova_results["n_per_condition"],
            "mean": anova_results["means"],
            "sem": anova_results["sems"],
        }
    )
    stats_df.to_csv(stats_path, index=False)

    # Save ANOVA summary
    anova_path = output_dir / f"{stem}_anova.txt"
    with open(anova_path, "w") as f:
        f.write(f"One-way ANOVA Results for {fig_name}\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"F-statistic: {anova_results['f_statistic']:.4f}\n")
        f.write(f"p-value: {anova_results['p_value']:.6f}\n\n")
        f.write("Group Statistics (using replicate means):\n")
        f.write("-" * 50 + "\n")
        for i, cond in enumerate(anova_results["conditions"]):
            f.write(
                f"  {cond}: n={anova_results['n_per_condition'][i]}, "
                f"mean={anova_results['means'][i]:.4f}, "
                f"SEM={anova_results['sems'][i]:.4f}\n"
            )
    logger.info(f"Saved ANOVA results to {anova_path}")

    return anova_results


def process_all_figures(
    figures_dir: Path,
    output_dir: Path | None = None,
    pattern: str = "*_normalized.csv",
) -> None:
    """
    Process all figure CSV files in a directory.

    Args:
        figures_dir: Directory containing figure CSV files
        output_dir: Output directory (defaults to figures_dir/plots)
        pattern: Glob pattern for CSV files to process
    """
    figures_dir = Path(figures_dir)
    output_dir = Path(output_dir) if output_dir else figures_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(figures_dir.glob(pattern))
    if not csv_files:
        logger.warning(f"No files matching {pattern} found in {figures_dir}")
        return

    logger.info(f"Found {len(csv_files)} files to process")

    all_results = []
    for csv_path in csv_files:
        # Determine ylabel based on metric type
        if "edge_spot" in csv_path.name:
            ylabel = "Edge Spot Fraction (normalized)"
        elif "gini" in csv_path.name:
            ylabel = "Gini Index (normalized)"
        else:
            ylabel = "Value (normalized)"

        results = process_figure_file(csv_path, output_dir, ylabel=ylabel)
        if results:
            results["file"] = csv_path.name
            all_results.append(results)

    # Save summary of all ANOVA results
    if all_results:
        summary_path = output_dir / "anova_summary.csv"
        summary_df = pd.DataFrame(
            [
                {
                    "file": r["file"],
                    "f_statistic": r["f_statistic"],
                    "p_value": r["p_value"],
                    "n_conditions": len(r["conditions"]),
                }
                for r in all_results
            ]
        )
        summary_df.to_csv(summary_path, index=False)
        logger.info(f"Saved ANOVA summary to {summary_path}")


def main():
    """CLI entry point for figure plotting."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Generate publication-ready figures with statistical analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all normalized figure CSVs
  edge-spot-plot -i results/figures/

  # Process specific file
  edge-spot-plot -i results/figures/Fig1D_edge_spot_normalized.csv

  # Output to specific directory
  edge-spot-plot -i results/figures/ -o plots/
        """,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input directory with figure CSVs or single CSV file",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="Output directory for plots (defaults to input/plots)",
    )

    parser.add_argument(
        "--pattern",
        "-p",
        default="*_normalized.csv",
        help="Glob pattern for CSV files (default: *_normalized.csv)",
    )

    args = parser.parse_args()

    input_path = args.input

    if input_path.is_file():
        # Process single file
        output_dir = args.output if args.output else input_path.parent / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        process_figure_file(input_path, output_dir)
    else:
        # Process directory
        process_all_figures(input_path, args.output, args.pattern)


if __name__ == "__main__":
    main()
