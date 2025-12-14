"""
I/O utilities for the CellProfiler pipeline port.

This module handles:
- File discovery and pattern matching
- Filename parsing and metadata extraction
- CSV export in CellProfiler-compatible format
- Image loading and normalization

Author: Generated from CellProfiler pipeline 231115_combined_pipeline_new_nomito_fixed_moments
Pipeline Version: v6
"""

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage import io

logger = logging.getLogger(__name__)


@dataclass
class ImagePair:
    """Represents a matched pair of Hoechst and MIRO160mer images."""

    hoechst_path: Path
    miro_path: Path
    well_number: str
    sequence: str
    xy_position: str
    date: str

    @property
    def image_id(self) -> str:
        """Unique identifier for this image pair."""
        return f"{self.date}_{self.well_number}_{self.sequence}_{self.xy_position}"

    def to_metadata_dict(self, image_number: int) -> dict[str, Any]:
        """
        Convert to metadata dictionary for CSV export.

        Args:
            image_number: Sequential image number

        Returns:
            Dictionary with CellProfiler-style metadata
        """
        return {
            'ImageNumber': image_number,
            'FileName_Hoechst': self.hoechst_path.name,
            'FileName_MIRO160mer': self.miro_path.name,
            'PathName_Hoechst': str(self.hoechst_path.parent),
            'PathName_MIRO160mer': str(self.miro_path.parent),
            'Metadata_Date': self.date,
            'Metadata_Well': self.well_number,
            'Metadata_Sequence': self.sequence,
            'Metadata_XY': self.xy_position,
        }


class FileDiscovery:
    """Discover and match image files from CellProfiler-compatible directory structure."""

    @staticmethod
    def find_image_pairs(
        input_dir: Path,
        date: str | None = None
    ) -> list[ImagePair]:
        """
        Find all Hoechst/MIRO160mer image pairs in a directory.

        Expected directory structure:
            input_dir/
                {date}/
                    {wellnumber}_{channel}_{sequence}_{xy}.tif

        Where:
            - channel is "405" for Hoechst or "488" for MIRO160mer
            - Files are matched by (wellnumber, sequence, xy)

        Args:
            input_dir: Root directory containing date subdirectories
            date: Optional specific date to process (if None, process all dates)

        Returns:
            List of ImagePair objects
        """
        input_path = Path(input_dir)

        if not input_path.exists():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")

        # Find all date directories
        if date is not None:
            date_dirs = [input_path / date]
            if not date_dirs[0].exists():
                raise FileNotFoundError(f"Date directory not found: {date_dirs[0]}")
        else:
            # Find all subdirectories that look like dates
            date_dirs = [d for d in input_path.iterdir() if d.is_dir()]

        image_pairs = []

        for date_dir in date_dirs:
            date_name = date_dir.name

            # Find all TIF files
            tif_files = list(date_dir.glob("*.tif")) + list(date_dir.glob("*.TIF"))

            # Separate by channel
            hoechst_files = [f for f in tif_files if "_405" in f.stem]
            miro_files = [f for f in tif_files if "_488" in f.stem]

            # Match pairs
            pairs = FileDiscovery._match_image_pairs(
                hoechst_files, miro_files, date_name
            )
            image_pairs.extend(pairs)

        return image_pairs

    @staticmethod
    def _match_image_pairs(
        hoechst_files: list[Path],
        miro_files: list[Path],
        date: str
    ) -> list[ImagePair]:
        """
        Match Hoechst and MIRO160mer files by their base names.

        Args:
            hoechst_files: List of Hoechst (405nm) file paths
            miro_files: List of MIRO160mer (488nm) file paths
            date: Date identifier

        Returns:
            List of matched ImagePair objects
        """
        pairs = []

        # Parse all Hoechst files
        hoechst_dict = {}
        for h_file in hoechst_files:
            parsed = FileDiscovery._parse_filename(h_file.stem)
            if parsed:
                key = (parsed['well'], parsed['sequence'], parsed['xy'])
                hoechst_dict[key] = (h_file, parsed)

        # Match with MIRO files
        for m_file in miro_files:
            parsed = FileDiscovery._parse_filename(m_file.stem)
            if parsed:
                key = (parsed['well'], parsed['sequence'], parsed['xy'])

                if key in hoechst_dict:
                    h_file, h_parsed = hoechst_dict[key]

                    pair = ImagePair(
                        hoechst_path=h_file,
                        miro_path=m_file,
                        well_number=parsed['well'],
                        sequence=parsed['sequence'],
                        xy_position=parsed['xy'],
                        date=date
                    )
                    pairs.append(pair)

        return pairs

    @staticmethod
    def _parse_filename(filename: str) -> dict[str, str] | None:
        """
        Parse filename to extract metadata.

        Expected format: Well{WELL}_Channel{...}_Seq{SEQ}-MaxIP_XY{XY}_{CHANNEL}
        Example: "WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_405"

        Args:
            filename: Filename without extension

        Returns:
            Dictionary with parsed components or None if parsing fails
        """
        # Pattern: [Plate000_]WellXXX_ChannelYYY_SeqZZZZ-MaxIP_XYAA_CCC
        # Optional Plate prefix to handle different naming conventions
        pattern = r"^(?:Plate\d+_)?Well([A-Z]\d+)_Channel[\d,]+_Seq(\d+)-MaxIP_XY(\d+)_(\d+)$"

        match = re.match(pattern, filename)
        if match:
            return {
                'well': match.group(1),
                'sequence': f"Seq{match.group(2)}",
                'xy': f"XY{match.group(3)}",
                'channel': match.group(4),
            }

        return None


class ImageLoader:
    """Load and normalize microscopy images."""

    @staticmethod
    def load_image(file_path: Path) -> np.ndarray:
        """
        Load a microscopy image from disk.

        Images are scaled to 0-1 range based on dtype (16-bit -> /65535).
        No min-max stretching is applied to preserve absolute intensity values.

        Args:
            file_path: Path to image file (TIF format)

        Returns:
            2D float32 array in 0-1 range
        """
        # Load image
        image = io.imread(str(file_path))

        # Convert to float
        if image.dtype != np.float64 and image.dtype != np.float32:
            # Rescale based on dtype
            if image.dtype == np.uint8:
                image = image.astype(np.float32) / 255.0
            elif image.dtype == np.uint16:
                image = image.astype(np.float32) / 65535.0
            else:
                image = image.astype(np.float32)

        # Ensure 2D
        if image.ndim > 2:
            # Take first channel or max projection
            image = image.max(axis=0) if image.shape[0] < image.shape[-1] else image[0]

        # Note: We don't apply min-max normalization here because:
        # 1. CellProfiler uses raw 0-65535 -> 0-1 scaling (dtype-based only)
        # 2. Min-max stretching destroys absolute intensity values needed for measurements
        # The dtype-based scaling above already converts to 0-1 range

        return image

    @staticmethod
    def load_image_pair(
        image_pair: ImagePair
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Load a matched pair of Hoechst and MIRO160mer images.

        Args:
            image_pair: ImagePair object with file paths

        Returns:
            Tuple of (hoechst_image, miro_image) as 2D float arrays
        """
        hoechst = ImageLoader.load_image(image_pair.hoechst_path)
        miro = ImageLoader.load_image(image_pair.miro_path)

        # Ensure same dimensions
        if hoechst.shape != miro.shape:
            raise ValueError(
                f"Image dimensions don't match: "
                f"Hoechst {hoechst.shape} vs MIRO {miro.shape}"
            )

        return hoechst, miro


class CSVExporter:
    """Export measurements to CellProfiler-compatible CSV files."""

    @staticmethod
    def export_measurements(
        measurements_list: list[dict[str, pd.DataFrame]],
        output_dir: Path,
        date: str | None = None
    ) -> None:
        """
        Export measurements to CSV files matching CellProfiler format.

        Args:
            measurements_list: List of measurement dictionaries (one per image)
                               Each dict has keys: 'Nuclei', 'Expand_Nuclei',
                               'Perinuclear_region', 'edge_spots', 'Image'
            output_dir: Output directory for CSV files
            date: Optional date identifier for subdirectory
        """
        output_path = Path(output_dir)
        if date:
            output_path = output_path / date
        output_path.mkdir(parents=True, exist_ok=True)

        # Combine measurements from all images
        combined = {
            'Nuclei': [],
            'Expand_Nuclei': [],
            'Perinuclear_region': [],
            'edge_spots': [],
            'Image': [],
        }

        for meas_dict in measurements_list:
            for key in combined.keys():
                if key in meas_dict and not meas_dict[key].empty:
                    combined[key].append(meas_dict[key])

        # Concatenate and export each object type
        for object_name, df_list in combined.items():
            if not df_list:
                continue

            combined_df = pd.concat(df_list, ignore_index=True)

            # Sort columns to match CellProfiler output
            combined_df = CSVExporter._sort_columns(combined_df)

            # Export to CSV
            output_file = output_path / f"{object_name}.csv"
            combined_df.to_csv(output_file, index=False)
            logger.info(f"Exported {len(combined_df)} rows to {output_file}")

        # Create All_measurements.csv (combined Nuclei and edge_spots)
        CSVExporter._create_all_measurements(combined, output_path)

    @staticmethod
    def _sort_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        Sort columns to match CellProfiler output format.

        Order:
        1. ImageNumber, ObjectNumber
        2. Metadata columns (FileName, PathName, Metadata_*)
        3. Number_* columns
        4. Measurement columns (alphabetically)
        """
        if df.empty:
            return df

        cols = df.columns.tolist()

        # Define column order
        priority_cols = []

        if 'ImageNumber' in cols:
            priority_cols.append('ImageNumber')
            cols.remove('ImageNumber')

        if 'ObjectNumber' in cols:
            priority_cols.append('ObjectNumber')
            cols.remove('ObjectNumber')

        # Metadata columns
        metadata_cols = [c for c in cols if c.startswith(('FileName_', 'PathName_', 'Metadata_'))]
        for c in metadata_cols:
            cols.remove(c)
        metadata_cols.sort()

        # Number columns
        number_cols = [c for c in cols if c.startswith('Number_')]
        for c in number_cols:
            cols.remove(c)
        number_cols.sort()

        # Remaining measurement columns (alphabetically)
        cols.sort()

        # Combine
        ordered_cols = priority_cols + metadata_cols + number_cols + cols

        return df[ordered_cols]

    @staticmethod
    def _create_all_measurements(
        combined: dict[str, list[pd.DataFrame]],
        output_path: Path
    ) -> None:
        """
        Create All_measurements.csv combining Nuclei and edge_spots.

        Produces CellProfiler-compatible multi-level column structure:
        - Level 0: Object type (Image, Nuclei, edge_spots)
        - Level 1: Column name (FileName_Hoechst, Number_Object_Number, etc.)

        Args:
            combined: Dictionary of combined measurements
            output_path: Output directory
        """
        if not combined['Nuclei'] or not combined['Image']:
            return

        nuclei_df = pd.concat(combined['Nuclei'], ignore_index=True)
        image_df = pd.concat(combined['Image'], ignore_index=True)

        # Handle edge_spots (may be empty for some images)
        if combined['edge_spots']:
            edge_spots_df = pd.concat(combined['edge_spots'], ignore_index=True)
        else:
            edge_spots_df = pd.DataFrame()

        # Build rows for each image with multi-level column structure
        rows = []
        for _, img_row in image_df.iterrows():
            image_number = img_row['ImageNumber']

            # Count nuclei for this image
            nuclei_count = len(nuclei_df[nuclei_df['ImageNumber'] == image_number])

            # Count edge spots for this image
            if not edge_spots_df.empty:
                edge_spot_count = len(
                    edge_spots_df[edge_spots_df['ImageNumber'] == image_number]
                )
            else:
                edge_spot_count = 0

            rows.append({
                ('Image', 'ImageNumber'): image_number,
                ('Image', 'FileName_Hoechst'): img_row.get('FileName_Hoechst', ''),
                ('Image', 'FileName_MIRO160mer'): img_row.get('FileName_MIRO160mer', ''),
                ('Nuclei', 'Number_Object_Number'): nuclei_count,
                ('edge_spots', 'Number_Object_Number'): edge_spot_count,
            })

        if rows:
            all_df = pd.DataFrame(rows)
            all_df.columns = pd.MultiIndex.from_tuples(all_df.columns)
            output_file = output_path / "All_measurements.csv"
            all_df.to_csv(output_file, index=False)
            logger.info(f"Exported {len(all_df)} rows to {output_file}")


def load_parameter_config(config_path: Path) -> dict[str, dict[str, Any]]:
    """
    Load per-date parameter configuration from JSON, YAML, or CSV file.

    Supported formats:

    JSON/YAML:
    {
        "default": {"otsu_correction_factor": 0.45, ...},
        "20231115": {"otsu_correction_factor": 0.50, ...}
    }

    CSV (simple table format):
    date,otsu_correction_factor,diameter_min,...
    default,0.45,30,...
    241223,0.50,,...
    241217,0.40,25,...

    Empty cells in CSV inherit from 'default' row.

    Args:
        config_path: Path to configuration file (JSON, YAML, or CSV)

    Returns:
        Dictionary mapping date to parameter overrides
    """
    config_path = Path(config_path)

    if config_path.suffix == '.csv':
        return _load_csv_config(config_path)
    elif config_path.suffix == '.json':
        import json
        with open(config_path) as f:
            return json.load(f)
    elif config_path.suffix in ['.yml', '.yaml']:
        import yaml
        with open(config_path) as f:
            return yaml.safe_load(f)
    else:
        raise ValueError(f"Unsupported config format: {config_path.suffix}")


def _load_csv_config(config_path: Path) -> dict[str, dict[str, Any]]:
    """
    Load parameter configuration from CSV file.

    CSV format:
    date,otsu_correction_factor,diameter_min,diameter_max,...
    default,0.45,30,100,...
    241223,0.50,,,...
    241217,0.40,25,,...

    Empty cells inherit from the 'default' row.
    The 'date' column is required and used as the key.

    Args:
        config_path: Path to CSV file

    Returns:
        Dictionary mapping date to parameter overrides
    """
    import csv

    config: dict[str, dict[str, Any]] = {}

    with open(config_path, newline='') as f:
        reader = csv.DictReader(f)

        if 'date' not in reader.fieldnames:
            raise ValueError("CSV config must have a 'date' column")

        for row in reader:
            date = row.pop('date').strip()
            if not date:
                continue

            # Convert values to appropriate types
            params: dict[str, Any] = {}
            for key, value in row.items():
                value = value.strip() if value else ''
                if not value:
                    continue  # Skip empty values (will inherit from default)

                # Try to convert to appropriate type
                params[key] = _parse_config_value(value)

            config[date] = params

    return config


def _parse_config_value(value: str) -> int | float | bool | str:
    """Parse a config value string to appropriate Python type."""
    # Try int
    try:
        return int(value)
    except ValueError:
        pass

    # Try float
    try:
        return float(value)
    except ValueError:
        pass

    # Try bool
    if value.lower() in ('true', 'yes', '1'):
        return True
    if value.lower() in ('false', 'no', '0'):
        return False

    # Return as string
    return value


def get_parameters_for_date(
    config: dict[str, dict[str, Any]],
    date: str
) -> dict[str, Any]:
    """
    Get parameters for a specific date, falling back to defaults.

    Args:
        config: Full parameter configuration
        date: Date identifier

    Returns:
        Dictionary of parameters for this date
    """
    # Start with defaults
    params = config.get('default', {}).copy()

    # Override with date-specific params
    if date in config:
        params.update(config[date])

    return params
