"""
Figure-based aggregation module.

This module creates aggregated tables for each Figure sheet in the configuration
Excel file, combining data across dates and wells into condition-based tables.
"""

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class FigureAggregator:
    """Aggregate data across dates/wells for figure generation."""

    def __init__(
        self,
        results_dir: Path,
        config_path: Path,
        output_dir: Path | None = None,
    ):
        """
        Initialize the figure aggregator.

        Args:
            results_dir: Directory containing date-based result folders
            config_path: Path to Excel config file with Figure sheets
            output_dir: Output directory (defaults to results_dir/figures)
        """
        self.results_dir = Path(results_dir)
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir) if output_dir else self.results_dir / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load Excel file
        self.excel = pd.ExcelFile(config_path)
        self.figure_sheets = [
            s for s in self.excel.sheet_names if s.startswith("Fig") and s not in ("Figures",)
        ]
        logger.info(f"Found {len(self.figure_sheets)} figure sheets: {self.figure_sheets}")

    def process_all_figures(self) -> None:
        """Process all figure sheets."""
        for sheet_name in self.figure_sheets:
            try:
                self.process_figure(sheet_name)
            except Exception as e:
                logger.error(f"Failed to process {sheet_name}: {e}")

    def process_figure(self, sheet_name: str) -> None:
        """
        Process a single figure sheet.

        Args:
            sheet_name: Name of the Excel sheet (e.g., 'Fig1D')
        """
        logger.info(f"Processing figure: {sheet_name}")

        # Parse the figure sheet
        fig_df = pd.read_excel(self.config_path, sheet_name=sheet_name, header=0)

        # First column is dates, rest are conditions
        date_col = fig_df.columns[0]
        conditions = [c for c in fig_df.columns[1:] if not str(c).startswith("Unnamed")]

        logger.info(f"Conditions: {conditions}")

        # Build data for edge spot fraction and Gini
        edge_spot_data = {}
        gini_data = {}

        for _, row in fig_df.iterrows():
            date_val = row[date_col]
            if pd.isna(date_val):
                continue

            # Convert date to string (handle both int and float)
            date_str = str(int(date_val))

            for condition in conditions:
                well_val = row[condition]
                if pd.isna(well_val):
                    continue

                # Handle comma-separated wells (e.g., "B06,D06,B07")
                wells = [w.strip() for w in str(well_val).split(",")]

                for well in wells:
                    col_name = f"{self._sanitize_name(condition)}_{date_str}_{well}"

                    # Try to load edge spot fraction data
                    edge_spot_col = self._load_well_data(
                        date_str, well, "edge_spot_fraction_static.csv"
                    )
                    if edge_spot_col is not None:
                        edge_spot_data[col_name] = edge_spot_col

                    # Try to load Gini data
                    gini_col = self._load_well_data(
                        date_str, well, "GINI_Gini_MIRO160mer_fov_median_static.csv"
                    )
                    if gini_col is not None:
                        gini_data[col_name] = gini_col

        # Create DataFrames and save
        if edge_spot_data:
            self._save_figure_outputs(edge_spot_data, sheet_name, "edge_spot", conditions)
        else:
            logger.warning(f"No edge spot data found for {sheet_name}")

        if gini_data:
            self._save_figure_outputs(gini_data, sheet_name, "gini", conditions)
        else:
            logger.warning(f"No Gini data found for {sheet_name}")

    def _load_well_data(self, date_str: str, well: str, filename: str) -> pd.Series | None:
        """
        Load a single well's column from a static CSV file.

        Args:
            date_str: Date string (e.g., '230912')
            well: Well identifier (e.g., 'B02')
            filename: Name of the CSV file to load

        Returns:
            Series with XY as index, or None if not found
        """
        file_path = self.results_dir / date_str / filename
        if not file_path.exists():
            logger.debug(f"File not found: {file_path}")
            return None

        try:
            df = pd.read_csv(file_path)

            # Handle files with XY column (Gini median files)
            if "XY" in df.columns:
                df = df.set_index("XY")
                # Keep only well columns (e.g., B02, D02, E02)
                well_cols = [c for c in df.columns if len(c) <= 3 and c[0].isupper()]
                if well_cols:
                    df = df[well_cols]
            else:
                # Edge spot files have XY as first column (index)
                df = pd.read_csv(file_path, index_col=0)

            if well not in df.columns:
                logger.debug(f"Well {well} not found in {file_path}")
                return None

            result = df[well]
            result.index.name = "XY"
            return result
        except Exception as e:
            logger.warning(f"Error loading {file_path}: {e}")
            return None

    def _save_figure_outputs(
        self,
        data: dict[str, pd.Series],
        fig_name: str,
        metric: str,
        conditions: list[str],
    ) -> None:
        """
        Save raw and normalized figure outputs.

        Args:
            data: Dict mapping column names to data series
            fig_name: Figure name (e.g., 'Fig1D')
            metric: Metric name ('edge_spot' or 'gini')
            conditions: List of condition names for normalization
        """
        # Create DataFrame with all columns
        df = pd.DataFrame(data)

        # Sort columns by condition order
        sorted_cols = self._sort_columns_by_condition(df.columns.tolist(), conditions)
        df = df[sorted_cols]

        # Save raw data
        raw_path = self.output_dir / f"{fig_name}_{metric}_raw.csv"
        df.to_csv(raw_path)
        logger.info(f"Wrote raw {metric} data to {raw_path}")

        # Calculate normalized data (normalize per-date by that date's control mean)
        # Column format: {condition}_{date}_{well}
        first_condition = self._sanitize_name(conditions[0])
        date_pattern = re.compile(r"_(\d{6})_")

        # Extract unique dates from column names
        dates = set()
        for col in df.columns:
            match = date_pattern.search(col)
            if match:
                dates.add(match.group(1))

        if dates:
            normalized_df = df.copy()
            for date in dates:
                date_cols = [c for c in df.columns if f"_{date}_" in c]
                # Match exact condition: {condition}_{date}_{well}
                # Use regex to avoid partial matches (e.g., "T2" matching "T2_E242R")
                control_pattern = re.compile(
                    rf"^{re.escape(first_condition)}_\d{{6}}_[A-H]\d{{2}}$"
                )
                control_cols = [c for c in date_cols if control_pattern.match(c)]

                if control_cols:
                    control_mean = df[control_cols].mean().mean()
                    if control_mean != 0:
                        normalized_df[date_cols] = df[date_cols] / control_mean
                    else:
                        logger.warning(
                            f"Control mean is zero for {fig_name} date {date}, "
                            "skipping normalization for this date"
                        )
                else:
                    logger.warning(f"No control columns found for {fig_name} date {date}")

            norm_path = self.output_dir / f"{fig_name}_{metric}_normalized.csv"
            normalized_df.to_csv(norm_path)
            logger.info(f"Wrote normalized {metric} data to {norm_path}")
        else:
            logger.warning(f"No dates found in columns for {fig_name}")

    def _sort_columns_by_condition(self, columns: list[str], conditions: list[str]) -> list[str]:
        """Sort columns by condition order, then by date/well."""
        sanitized_conditions = [self._sanitize_name(c) for c in conditions]

        def sort_key(col: str) -> tuple:
            for i, cond in enumerate(sanitized_conditions):
                if col.startswith(f"{cond}_"):
                    return (i, col)
            return (len(conditions), col)

        return sorted(columns, key=sort_key)

    @staticmethod
    def _sanitize_name(name: str) -> str:
        """Sanitize condition name for use in column names."""
        return str(name).replace(" ", "_").replace("-", "_")


def main():
    """CLI entry point for figure aggregation."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Aggregate data across dates/wells for figure generation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all figures using default paths
  edge-spot-figure-aggregate -r results/ -c config.xlsx

  # Process to specific output directory
  edge-spot-figure-aggregate -r results/ -c config.xlsx -o figures/
        """,
    )

    parser.add_argument(
        "--results",
        "-r",
        type=Path,
        required=True,
        help="Directory containing date-based result folders",
    )

    parser.add_argument(
        "--config",
        "-c",
        type=Path,
        required=True,
        help="Excel config file with Figure sheets",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="Output directory (defaults to results/figures)",
    )

    args = parser.parse_args()

    aggregator = FigureAggregator(
        results_dir=args.results,
        config_path=args.config,
        output_dir=args.output,
    )
    aggregator.process_all_figures()


if __name__ == "__main__":
    main()
