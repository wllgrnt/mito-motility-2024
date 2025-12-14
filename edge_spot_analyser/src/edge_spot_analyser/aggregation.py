"""
Aggregation module for processing pipeline outputs.

This module transforms raw per-object CSV files into aggregated
plottable formats (edge spot fractions, mass displacement, CoV matrices).

Supports both static mode (single timepoint) and time-course mode.
"""

import argparse
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class AggregationConfig:
    """Configuration for aggregation processing."""

    # File names
    edge_spot_file: str = "All_measurements.csv"
    mass_displacement_file: str = "Expand_Nuclei.csv"
    cov_file: str = "Perinuclear_region.csv"

    # Processing options
    t_varies: bool = False
    plot: bool = False

    # Column mappings (support multiple possible column names for flexibility)
    mass_displacement_cols: dict[str, list[str]] = field(default_factory=lambda: {
        "mass_displacement_60mer": [
            "Intensity_MassDisplacement_MIRO160mer",
            "Intensity_MassDisplacement_MIRO160mer_rescaled",
        ],
    })

    # Extra columns to extract alongside CoV
    bonus_cols: list[str] = field(default_factory=lambda: [
        "GINI_Gini_MIRO160mer",
    ])


class Aggregator:
    """Process pipeline outputs into aggregated analysis files."""

    def __init__(self, config: AggregationConfig | None = None):
        self.config = config or AggregationConfig()

    def process_date(
        self,
        input_dir: Path,
        output_dir: Path | None = None,
    ) -> None:
        """
        Process a single date's pipeline output.

        Args:
            input_dir: Directory containing CSV files for this date
            output_dir: Directory for aggregated outputs (defaults to input_dir)
        """
        input_dir = Path(input_dir)
        output_dir = Path(output_dir) if output_dir else input_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Processing aggregation for {input_dir}")

        # Edge spot fraction
        edge_spot_file = input_dir / self.config.edge_spot_file
        if edge_spot_file.exists():
            self._process_edge_spots(edge_spot_file, output_dir)
        else:
            logger.warning(f"Edge spot file not found: {edge_spot_file}")

        # Mass displacement
        mass_file = input_dir / self.config.mass_displacement_file
        if mass_file.exists():
            self._process_mass_displacement(mass_file, output_dir)
        else:
            logger.warning(f"Mass displacement file not found: {mass_file}")

        # CoV and bonus columns
        cov_file = input_dir / self.config.cov_file
        if cov_file.exists():
            self._process_cov(cov_file, output_dir)
            self._process_bonus_cols(cov_file, output_dir)
        else:
            logger.warning(f"CoV file not found: {cov_file}")

    def _process_edge_spots(self, input_path: Path, output_dir: Path) -> None:
        """
        Generate edge_spot_count / nuclei_count for each field of view.

        Args:
            input_path: Path to All_measurements.csv
            output_dir: Output directory for aggregated files
        """
        logger.info(f"Processing edge spots from {input_path}")

        # Read with multi-level headers
        raw_df = pd.read_csv(input_path, header=[0, 1])
        processed_df = self._extract_edgespot_cols(raw_df)

        # Save intermediate raw data
        intermediate_path = output_dir / "edge_spot_fraction_raw.csv"
        processed_df.to_csv(intermediate_path, index=False)
        logger.info(f"Wrote edge spot intermediate to {intermediate_path}")

        if self.config.t_varies:
            self._generate_edge_spot_time_files(processed_df, output_dir)
        else:
            self._generate_edge_spot_static_files(processed_df, output_dir)

    def _extract_edgespot_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract well number, XY, T and aggregate to get edge spot fraction."""
        logger.info(f"Extracting edge spot columns, input shape {df.shape}")

        filename_column = "FileName_Hoechst"
        nuclei_count_column = "Number_Object_Number"
        edge_spot_count_column = "Number_Object_Number"

        output_df = df.copy()

        if self.config.t_varies:
            output_df["T"] = output_df["Image"][filename_column].apply(extract_timestamp)
        output_df["XY"] = output_df["Image"][filename_column].apply(extract_xy)
        output_df["WellNumber"] = output_df["Image"][filename_column].apply(
            extract_wellnumber
        )
        output_df["nuclei_count"] = output_df["Nuclei"][nuclei_count_column]
        output_df["edge_spot_count"] = output_df["edge_spots"][edge_spot_count_column]

        cols = ["WellNumber", "XY", "nuclei_count", "edge_spot_count"]
        if self.config.t_varies:
            cols.append("T")

        output_df = output_df[cols].reset_index(drop=True)

        # Drop multi-index level if present
        if isinstance(output_df.columns, pd.MultiIndex):
            output_df.columns = output_df.columns.droplevel(1)

        # Calculate edge spot fraction
        output_df["edge_spot_fraction"] = (
            output_df["edge_spot_count"] / output_df["nuclei_count"]
        )

        logger.info(f"Extracted edge spot data, output shape {output_df.shape}")
        return output_df

    def _generate_edge_spot_time_files(
        self, df: pd.DataFrame, output_dir: Path
    ) -> None:
        """Generate time-course edge spot files."""
        # File for each well number, T as columns, XY as rows
        for well_number in df.WellNumber.unique():
            well_df = df[df.WellNumber == well_number]
            output_path = output_dir / f"edge_spot_fraction_over_time_{well_number}.csv"
            self._save_pivot_table(
                data=well_df,
                values="edge_spot_fraction",
                index="XY",
                columns="T",
                output_path=output_path,
            )

        # Average over XYs, normalise and optionally plot
        average_over_time = (
            df.groupby(["WellNumber", "T"]).edge_spot_fraction.mean().reset_index()
        )
        pivot = pd.pivot_table(
            data=average_over_time,
            values="edge_spot_fraction",
            index="T",
            columns="WellNumber",
        )
        pivot = pivot.divide(pivot.iloc[0])
        output_path = output_dir / "edge_spot_fraction_mean_over_time_normalised.csv"
        logger.info(f"Writing normalised pivot table to {output_path}")
        pivot.to_csv(output_path)

        if self.config.plot:
            self._plot_time_series(pivot, "Edge spot fraction over time, normalised to T0")

    def _generate_edge_spot_static_files(
        self, df: pd.DataFrame, output_dir: Path
    ) -> None:
        """Generate static (single timepoint) edge spot files."""
        output_path = output_dir / "edge_spot_fraction_static.csv"
        self._save_pivot_table(
            data=df,
            values="edge_spot_fraction",
            index="XY",
            columns="WellNumber",
            output_path=output_path,
        )

    def _process_mass_displacement(self, input_path: Path, output_dir: Path) -> None:
        """Process mass displacement data."""
        logger.info(f"Processing mass displacement from {input_path}")

        raw_df = pd.read_csv(input_path)
        processed_df = self._extract_massdisplacement_cols(raw_df)

        for displacement_type in self.config.mass_displacement_cols:
            if displacement_type in processed_df.columns:
                self._generate_ragged_df(
                    processed_df,
                    data_column=displacement_type,
                    output_dir=output_dir,
                )

    def _extract_massdisplacement_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract well number, XY, T, and mass displacement columns."""
        logger.info(f"Extracting mass displacement columns, input shape {df.shape}")

        filename_column = "FileName_MIRO160mer"
        output_df = df.copy()

        if self.config.t_varies:
            output_df["T"] = output_df[filename_column].apply(extract_timestamp)
        output_df["XY"] = output_df[filename_column].apply(extract_xy)
        output_df["WellNumber"] = output_df[filename_column].apply(extract_wellnumber)

        cols = ["WellNumber", "XY"]
        if self.config.t_varies:
            cols.append("T")

        # Find mass displacement columns
        for name, possible_cols in self.config.mass_displacement_cols.items():
            for col in possible_cols:
                if col in df.columns:
                    output_df[name] = df[col]
                    cols.append(name)
                    break

        return output_df[cols].reset_index(drop=True)

    def _process_cov(self, input_path: Path, output_dir: Path) -> None:
        """Process coefficient of variation data."""
        logger.info(f"Processing CoV from {input_path}")

        raw_df = pd.read_csv(input_path)
        processed_df = self._extract_cov_cols(raw_df)

        self._generate_ragged_df(
            processed_df,
            data_column="CoV",
            output_dir=output_dir,
        )

    def _extract_cov_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract well number, XY, T, and calculate CoV."""
        logger.info(f"Extracting CoV columns, input shape {df.shape}")

        filename_column = "FileName_MIRO160mer"
        std_columns = [
            "Intensity_StdIntensity_MIRO160mer",
            "Intensity_StdIntensity_MIRO160mer_rescaled",
        ]
        mean_columns = [
            "Intensity_MeanIntensity_MIRO160mer",
            "Intensity_MeanIntensity_MIRO160mer_rescaled",
        ]

        output_df = df.copy()

        if self.config.t_varies:
            output_df["T"] = output_df[filename_column].apply(extract_timestamp)
        output_df["XY"] = output_df[filename_column].apply(extract_xy)
        output_df["WellNumber"] = output_df[filename_column].apply(extract_wellnumber)

        # Calculate CoV
        for std_col, mean_col in zip(std_columns, mean_columns, strict=True):
            if std_col in df.columns and mean_col in df.columns:
                output_df["CoV"] = df[std_col] / df[mean_col]
                break
        else:
            raise ValueError("Could not find columns for CoV calculation")

        cols = ["WellNumber", "XY", "CoV"]
        if self.config.t_varies:
            cols.append("T")

        return output_df[cols].reset_index(drop=True)

    def _process_bonus_cols(self, input_path: Path, output_dir: Path) -> None:
        """Process additional dispersion measure columns."""
        if not self.config.bonus_cols:
            return

        logger.info(f"Processing bonus columns from {input_path}")

        raw_df = pd.read_csv(input_path)

        filename_column = "FileName_MIRO160mer"
        output_df = raw_df.copy()

        if self.config.t_varies:
            output_df["T"] = output_df[filename_column].apply(extract_timestamp)
        output_df["XY"] = output_df[filename_column].apply(extract_xy)
        output_df["WellNumber"] = output_df[filename_column].apply(extract_wellnumber)

        for bonus_col in self.config.bonus_cols:
            if bonus_col in raw_df.columns:
                cols = ["WellNumber", "XY", bonus_col]
                if self.config.t_varies:
                    cols.append("T")
                bonus_df = output_df[cols].reset_index(drop=True)
                self._generate_ragged_df(
                    bonus_df,
                    data_column=bonus_col,
                    output_dir=output_dir,
                )

    def _generate_ragged_df(
        self,
        input_df: pd.DataFrame,
        data_column: str,
        output_dir: Path,
    ) -> None:
        """Generate ragged dataframes for per-cell measurements."""
        if self.config.t_varies:
            self._generate_ragged_time_files(input_df, data_column, output_dir)
        else:
            self._generate_ragged_static_files(input_df, data_column, output_dir)

    def _generate_ragged_time_files(
        self,
        df: pd.DataFrame,
        data_column: str,
        output_dir: Path,
    ) -> None:
        """Generate time-course ragged dataframes."""
        subdf = df[["WellNumber", "XY", "T", data_column]]

        for well_number in subdf.WellNumber.unique():
            well_df = subdf[subdf.WellNumber == well_number]

            # Complex multi-level structure with T and XY as columns
            output_path = output_dir / f"{data_column}_time_and_xy_{well_number}.csv"
            temp_df = (
                well_df.drop(columns=["WellNumber"]).set_index(["T", "XY"]).squeeze()
            )
            temp_df = temp_df.reset_index(name="value")
            temp_df["index"] = temp_df.groupby(["T", "XY"]).cumcount()
            temp_df = temp_df.set_index(["T", "XY", "index"])["value"].unstack(
                level=[0, 1]
            )
            logger.info(f"Writing {data_column} table with shape {temp_df.shape}")
            temp_df.to_csv(output_path)

            # Stacked variant: column for each T
            stacked_df = well_df.drop(columns=["WellNumber", "XY"])
            stacked_df["count"] = stacked_df.groupby("T").cumcount()
            stacked_df = stacked_df.pivot(index="count", columns="T", values=data_column)
            stacked_path = output_dir / f"{data_column}_stacked_{well_number}.csv"
            logger.info(f"Writing stacked {data_column} table with shape {stacked_df.shape}")
            stacked_df.to_csv(stacked_path)

        # Average over Well, T, save and optionally plot
        average_over_time = (
            subdf.groupby(["WellNumber", "T"])[data_column].mean().reset_index()
        )
        pivot = pd.pivot_table(
            data=average_over_time,
            values=data_column,
            index="T",
            columns="WellNumber",
        )
        pivot = pivot.divide(pivot.iloc[0])
        output_path = output_dir / f"{data_column}_mean_over_time_normalised.csv"
        logger.info(f"Writing normalised pivot table to {output_path}")
        pivot.to_csv(output_path)

        if self.config.plot:
            self._plot_time_series(pivot, f"Mean {data_column} over time, normalised to T0")

    def _generate_ragged_static_files(
        self,
        df: pd.DataFrame,
        data_column: str,
        output_dir: Path,
    ) -> None:
        """Generate static ragged dataframes."""
        # Table with WellNumber and XY as columns
        output_path = output_dir / f"{data_column}_static.csv"
        subdf = df[["WellNumber", "XY", data_column]].copy()
        subdf["index"] = subdf.groupby(["WellNumber", "XY"]).cumcount()
        subdf = subdf.pivot(index="index", columns=["WellNumber", "XY"], values=data_column)
        logger.info(f"Writing static {data_column} table with shape {subdf.shape}")
        subdf.to_csv(output_path)

        # Table with median per FoV
        median_df = df.groupby(["WellNumber", "XY"])[[data_column]].median().unstack().T
        median_path = output_dir / f"{data_column}_fov_median_static.csv"
        median_df.to_csv(median_path)

        # Stacked variant: column for each WellNumber
        stacked_df = df[["WellNumber", data_column]].copy()
        stacked_df["count"] = stacked_df.groupby("WellNumber").cumcount()
        stacked_df = stacked_df.pivot(index="count", columns="WellNumber", values=data_column)
        stacked_path = output_dir / f"{data_column}_stacked_static.csv"
        logger.info(f"Writing stacked {data_column} table with shape {stacked_df.shape}")
        stacked_df.to_csv(stacked_path)

    def _save_pivot_table(
        self,
        data: pd.DataFrame,
        values: str,
        index: str,
        columns: str,
        output_path: Path,
    ) -> None:
        """Create and save a pivot table."""
        pivot = pd.pivot_table(data=data, values=values, index=index, columns=columns)
        logger.info(f"Writing pivot table to {output_path}")
        pivot.to_csv(output_path)

    def _plot_time_series(self, df: pd.DataFrame, title: str) -> None:
        """Plot time series data if plotting is enabled."""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            sns.set_style("whitegrid")
            df.plot(title=title)
            plt.show()
        except ImportError:
            logger.warning(
                "matplotlib/seaborn not installed. Install with: pip install edge-spot-analyser[plot]"
            )


def extract_xy(input_string: str) -> int:
    """Extract XY position from filename."""
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        # Fallback pattern for different naming conventions
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    """Extract timestamp (T) from filename."""
    match = re.search(r"(?<=_T)\d+", input_string)
    if match:
        return int(match.group(0))
    return 0


def extract_wellnumber(input_string: str) -> str:
    """Extract well number from filename."""
    match = re.search(r"(?<=Well)[A-Z]\d+", input_string)
    if match:
        return match.group(0)
    raise ValueError(f"Could not extract well number from: {input_string}")


def main():
    """CLI entry point for standalone aggregation."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Aggregate pipeline outputs into analysis files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Aggregate outputs in place
  edge-spot-aggregate -i results/241223/

  # Aggregate to different output directory
  edge-spot-aggregate -i results/241223/ -o analysis/241223/

  # Enable time-course mode
  edge-spot-aggregate -i results/241223/ --t-varies

  # Generate plots
  edge-spot-aggregate -i results/241223/ --plot
        """,
    )

    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input directory containing CSV files",
    )

    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output directory for aggregated files (defaults to input directory)",
    )

    parser.add_argument(
        "--t-varies",
        action="store_true",
        help="Enable time-course mode (multiple timepoints)",
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots (requires matplotlib and seaborn)",
    )

    args = parser.parse_args()

    config = AggregationConfig(t_varies=args.t_varies, plot=args.plot)
    aggregator = Aggregator(config)
    aggregator.process_date(args.input, args.output)


if __name__ == "__main__":
    main()
