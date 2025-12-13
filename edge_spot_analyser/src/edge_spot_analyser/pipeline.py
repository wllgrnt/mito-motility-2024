"""
Main pipeline script for CellProfiler pipeline port.

This script orchestrates the complete image analysis pipeline:
1. Load Hoechst and MIRO160mer image pairs
2. Segment nuclei (3-class Otsu)
3. Create perinuclear regions
4. Detect edge spots (Robust Background)
5. Measure all properties
6. Export to CellProfiler-compatible CSVs

Usage:
    python pipeline.py --input inputs/ --output results/ --config params.json
    python pipeline.py --input inputs/ --output results/ --date 20231115

Author: Generated from CellProfiler pipeline 231115_combined_pipeline_new_nomito_fixed_moments
Pipeline Version: v6
"""

import argparse
import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from edge_spot_analyser.io_utils import (
    CSVExporter,
    FileDiscovery,
    ImageLoader,
    ImagePair,
    get_parameters_for_date,
    load_parameter_config,
)
from edge_spot_analyser.measurements import combine_measurements_for_export
from edge_spot_analyser.segmentation import (
    EdgeSpotParams,
    NucleiSegmentationParams,
    PerinuclearRegionParams,
    create_perinuclear_regions,
    detect_edge_spots,
    filter_edge_spots_by_edge_intensity,
    mask_peripheral_regions,
    segment_nuclei,
)

logger = logging.getLogger(__name__)


def _process_single_image(args: tuple) -> dict[str, Any]:
    """
    Process a single image pair. Module-level function for multiprocessing.

    Args:
        args: Tuple of (image_pair, image_number, nuclei_params, edge_spot_params, perinuclear_params)

    Returns:
        Dictionary with 'image_id', 'measurements' or 'error'
    """
    image_pair, image_number, nuclei_params, edge_spot_params, perinuclear_params = args

    try:
        # Load images
        hoechst, miro = ImageLoader.load_image_pair(image_pair)

        # Segment nuclei
        nuclei_labels = segment_nuclei(hoechst, nuclei_params)
        n_nuclei = nuclei_labels.max()

        if n_nuclei == 0:
            metadata = image_pair.to_metadata_dict(image_number)
            return {
                'image_id': image_pair.image_id,
                'measurements': {
                    'Nuclei': pd.DataFrame(),
                    'Expand_Nuclei': pd.DataFrame(),
                    'Perinuclear_region': pd.DataFrame(),
                    'edge_spots': pd.DataFrame(),
                    'Image': pd.DataFrame([metadata]),
                }
            }

        # Create perinuclear regions
        expand_nuclei_10px, expand_nuclei_15px, perinuclear_ring = \
            create_perinuclear_regions(nuclei_labels, perinuclear_params)

        # Mask peripheral regions
        masked_miro = mask_peripheral_regions(miro, expand_nuclei_15px)

        # Detect edge spots
        edge_spots_labels = detect_edge_spots(masked_miro, edge_spot_params)
        n_edge_spots = edge_spots_labels.max()

        # Filter edge spots
        if n_edge_spots > 0:
            edge_spots_labels = filter_edge_spots_by_edge_intensity(
                edge_spots_labels, miro, min_edge_intensity=0.0, max_edge_intensity=1.0
            )

        # Measure properties
        metadata = image_pair.to_metadata_dict(image_number)
        measurements = combine_measurements_for_export(
            nuclei_labels=nuclei_labels,
            expand_nuclei_labels=expand_nuclei_10px,
            perinuclear_region_labels=perinuclear_ring,
            edge_spots_labels=edge_spots_labels,
            hoechst_image=hoechst,
            miro_image=miro,
            metadata=metadata
        )

        return {
            'image_id': image_pair.image_id,
            'measurements': measurements,
            'n_nuclei': n_nuclei,
            'n_edge_spots': edge_spots_labels.max()
        }

    except Exception as e:
        return {
            'image_id': image_pair.image_id,
            'error': str(e)
        }


class Pipeline:
    """Complete image analysis pipeline."""

    def __init__(
        self,
        nuclei_params: Optional[NucleiSegmentationParams] = None,
        edge_spot_params: Optional[EdgeSpotParams] = None,
        perinuclear_params: Optional[PerinuclearRegionParams] = None
    ):
        """
        Initialize pipeline with parameters.

        Args:
            nuclei_params: Nuclei segmentation parameters (uses defaults if None)
            edge_spot_params: Edge spot detection parameters (uses defaults if None)
            perinuclear_params: Perinuclear region parameters (uses defaults if None)
        """
        self.nuclei_params = nuclei_params or NucleiSegmentationParams()
        self.edge_spot_params = edge_spot_params or EdgeSpotParams()
        self.perinuclear_params = perinuclear_params or PerinuclearRegionParams()

    def process_image_pair(
        self,
        image_pair: ImagePair,
        image_number: int
    ) -> dict[str, pd.DataFrame]:
        """
        Process a single image pair through the complete pipeline.

        Args:
            image_pair: ImagePair object with file paths
            image_number: Sequential image number

        Returns:
            Dictionary of measurements for all object types
        """
        logger.info(f"Processing {image_pair.image_id}")

        # Load images
        hoechst, miro = ImageLoader.load_image_pair(image_pair)

        # Step 1: Segment nuclei (Module 7)
        logger.debug("  Segmenting nuclei...")
        nuclei_labels = segment_nuclei(hoechst, self.nuclei_params)
        n_nuclei = nuclei_labels.max()
        logger.debug(f"  Found {n_nuclei} nuclei")

        if n_nuclei == 0:
            logger.warning(f"  No nuclei found in {image_pair.image_id}, skipping")
            return self._empty_measurements(image_pair, image_number)

        # Step 2: Create perinuclear regions (Modules 8-10, 14)
        logger.debug("  Creating perinuclear regions...")
        expand_nuclei_10px, expand_nuclei_15px, perinuclear_ring = \
            create_perinuclear_regions(nuclei_labels, self.perinuclear_params)

        # Step 3: Mask peripheral regions (Module 10)
        logger.debug("  Masking peripheral regions...")
        masked_miro = mask_peripheral_regions(miro, expand_nuclei_15px)

        # Step 4: Detect edge spots (Module 11)
        logger.debug("  Detecting edge spots...")
        edge_spots_labels = detect_edge_spots(masked_miro, self.edge_spot_params)
        n_edge_spots = edge_spots_labels.max()
        logger.debug(f"  Found {n_edge_spots} edge spots")

        # Step 5: Filter edge spots by edge intensity (Module 13)
        if n_edge_spots > 0:
            logger.debug("  Filtering edge spots by edge intensity...")
            edge_spots_labels = filter_edge_spots_by_edge_intensity(
                edge_spots_labels,
                miro,
                min_edge_intensity=0.0,
                max_edge_intensity=1.0
            )
            n_edge_spots_filtered = edge_spots_labels.max()
            logger.debug(f"  {n_edge_spots_filtered} spots after filtering")

        # Step 6: Measure all properties (Modules 12-20)
        logger.debug("  Measuring properties...")
        metadata = image_pair.to_metadata_dict(image_number)

        measurements = combine_measurements_for_export(
            nuclei_labels=nuclei_labels,
            expand_nuclei_labels=expand_nuclei_10px,
            perinuclear_region_labels=perinuclear_ring,
            edge_spots_labels=edge_spots_labels,
            hoechst_image=hoechst,
            miro_image=miro,
            metadata=metadata
        )

        logger.info(
            f"  Complete: {n_nuclei} nuclei, {edge_spots_labels.max()} edge spots"
        )

        return measurements

    def _empty_measurements(
        self,
        image_pair: ImagePair,
        image_number: int
    ) -> dict[str, pd.DataFrame]:
        """
        Create empty measurement dictionary for images with no nuclei.

        Args:
            image_pair: ImagePair object
            image_number: Sequential image number

        Returns:
            Dictionary with empty DataFrames
        """
        metadata = image_pair.to_metadata_dict(image_number)

        return {
            'Nuclei': pd.DataFrame(),
            'Expand_Nuclei': pd.DataFrame(),
            'Perinuclear_region': pd.DataFrame(),
            'edge_spots': pd.DataFrame(),
            'Image': pd.DataFrame([metadata]),
        }

    def process_date_batch(
        self,
        input_dir: Path,
        output_dir: Path,
        date: str,
        n_workers: int = 1
    ) -> None:
        """
        Process all images for a specific date.

        Args:
            input_dir: Input directory containing date subdirectories
            output_dir: Output directory for CSV files
            date: Date identifier
            n_workers: Number of parallel workers (1 = sequential)
        """
        logger.info(f"Processing date: {date}")

        # Find all image pairs for this date
        image_pairs = FileDiscovery.find_image_pairs(input_dir, date=date)

        if not image_pairs:
            logger.warning(f"No image pairs found for date {date}")
            return

        logger.info(f"Found {len(image_pairs)} image pairs")

        # Prepare args for parallel processing
        args_list = [
            (pair, i, self.nuclei_params, self.edge_spot_params, self.perinuclear_params)
            for i, pair in enumerate(image_pairs, start=1)
        ]

        all_measurements = []

        if n_workers == 1:
            # Sequential processing
            for args in args_list:
                result = _process_single_image(args)
                if 'error' in result:
                    logger.error(f"Error processing {result['image_id']}: {result['error']}")
                else:
                    logger.info(f"  {result['image_id']}: {result.get('n_nuclei', 0)} nuclei, {result.get('n_edge_spots', 0)} edge spots")
                    all_measurements.append(result['measurements'])
        else:
            # Parallel processing
            logger.info(f"Using {n_workers} parallel workers")
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {executor.submit(_process_single_image, args): args[0].image_id for args in args_list}

                for future in as_completed(futures):
                    image_id = futures[future]
                    try:
                        result = future.result()
                        if 'error' in result:
                            logger.error(f"Error processing {result['image_id']}: {result['error']}")
                        else:
                            logger.info(f"  {result['image_id']}: {result.get('n_nuclei', 0)} nuclei, {result.get('n_edge_spots', 0)} edge spots")
                            all_measurements.append(result['measurements'])
                    except Exception as e:
                        logger.error(f"Error processing {image_id}: {e}")

        # Export to CSV
        if all_measurements:
            logger.info(f"Exporting measurements to {output_dir}")
            CSVExporter.export_measurements(all_measurements, output_dir, date=date)
        else:
            logger.warning("No measurements to export")

    def process_all_dates(
        self,
        input_dir: Path,
        output_dir: Path,
        config: Optional[dict[str, dict[str, Any]]] = None,
        n_workers: int = 1
    ) -> None:
        """
        Process all dates in the input directory.

        Args:
            input_dir: Input directory containing date subdirectories
            output_dir: Output directory for CSV files
            config: Optional per-date parameter configuration
            n_workers: Number of parallel workers
        """
        input_path = Path(input_dir)

        # Find all date directories
        date_dirs = [d for d in input_path.iterdir() if d.is_dir()]

        if not date_dirs:
            logger.error(f"No date directories found in {input_dir}")
            return

        logger.info(f"Found {len(date_dirs)} date directories")

        for date_dir in date_dirs:
            date = date_dir.name

            # Update parameters for this date if config provided
            if config:
                date_params = get_parameters_for_date(config, date)
                self._update_parameters(date_params)
                logger.info(f"Using parameters for date {date}: {date_params}")

            # Process this date
            try:
                self.process_date_batch(input_dir, output_dir, date, n_workers=n_workers)
            except Exception as e:
                logger.error(f"Error processing date {date}: {e}")
                logger.exception(e)
                continue

    def _update_parameters(self, params: dict[str, Any]) -> None:
        """
        Update pipeline parameters from configuration.

        Args:
            params: Dictionary of parameter overrides
        """
        # Update nuclei segmentation params
        if 'otsu_correction_factor' in params:
            self.nuclei_params.otsu_correction_factor = params['otsu_correction_factor']
        if 'diameter_min' in params:
            self.nuclei_params.diameter_min = params['diameter_min']
        if 'diameter_max' in params:
            self.nuclei_params.diameter_max = params['diameter_max']
        if 'min_distance' in params:
            self.nuclei_params.min_distance = params['min_distance']

        # Update edge spot params
        if 'edge_spot_diameter_min' in params:
            self.edge_spot_params.diameter_min = params['edge_spot_diameter_min']
        if 'edge_spot_diameter_max' in params:
            self.edge_spot_params.diameter_max = params['edge_spot_diameter_max']
        if 'edge_spot_correction_factor' in params:
            self.edge_spot_params.correction_factor = params['edge_spot_correction_factor']

        # Update perinuclear region params
        if 'inner_expansion' in params:
            self.perinuclear_params.inner_expansion = params['inner_expansion']
        if 'outer_expansion' in params:
            self.perinuclear_params.outer_expansion = params['outer_expansion']


def main():
    """Main entry point for CLI."""
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(
        description='Run CellProfiler pipeline in Python',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all dates with default parameters
  python pipeline.py --input inputs/ --output results/

  # Process specific date
  python pipeline.py --input inputs/ --output results/ --date 20231115

  # Use custom parameter configuration
  python pipeline.py --input inputs/ --output results/ --config params.json

  # Verbose logging
  python pipeline.py --input inputs/ --output results/ --verbose
        """
    )

    parser.add_argument(
        '--input', '-i',
        type=Path,
        required=True,
        help='Input directory containing date subdirectories with TIF files'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        required=True,
        help='Output directory for CSV files'
    )

    parser.add_argument(
        '--date', '-d',
        type=str,
        help='Process only this specific date (optional)'
    )

    parser.add_argument(
        '--config', '-c',
        type=Path,
        help='Path to parameter configuration file (JSON or YAML)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )

    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=1,
        help='Number of parallel workers (default: 1, use -1 for all CPUs)'
    )

    args = parser.parse_args()

    # Handle workers
    n_workers = args.workers
    if n_workers == -1:
        n_workers = os.cpu_count() or 1
    elif n_workers < 1:
        n_workers = 1

    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load configuration if provided
    config = None
    if args.config:
        logger.info(f"Loading configuration from {args.config}")
        config = load_parameter_config(args.config)

    # Create pipeline
    pipeline = Pipeline()

    # Process images
    start_time = time.time()

    if args.date:
        # Process single date
        pipeline.process_date_batch(args.input, args.output, args.date, n_workers=n_workers)
    else:
        # Process all dates
        pipeline.process_all_dates(args.input, args.output, config, n_workers=n_workers)

    elapsed_time = time.time() - start_time
    logger.info(f"Pipeline complete in {elapsed_time:.2f} seconds")


if __name__ == '__main__':
    main()
