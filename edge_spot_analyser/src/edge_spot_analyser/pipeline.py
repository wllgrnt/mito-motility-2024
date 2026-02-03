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
from typing import Any

import pandas as pd

from edge_spot_analyser.io_utils import (
    CSVExporter,
    FileDiscovery,
    ImageLoader,
    ImagePair,
    extract_figure_inclusions,
    get_parameters_for_date,
    load_parameter_config,
    should_exclude_image,
    should_include_image,
)
from edge_spot_analyser.measurements import combine_measurements_for_export
from edge_spot_analyser.segmentation import (
    EdgeSpotParams,
    NucleiSegmentationParams,
    PerinuclearRegionParams,
    SmoothParams,
    create_perinuclear_regions,
    detect_edge_spots,
    filter_edge_spots_by_edge_intensity,
    filter_edge_spots_by_nuclei_proximity,
    mask_peripheral_regions,
    segment_nuclei_keep_border,
    smooth_image,
)

logger = logging.getLogger(__name__)


def _process_single_image(args: tuple) -> dict[str, Any]:
    """
    Process a single image pair. Module-level function for multiprocessing.

    Args:
        args: Tuple of (image_pair, image_number, smooth_params, nuclei_params, edge_spot_params, perinuclear_params)

    Returns:
        Dictionary with 'image_id', 'measurements' or 'error'
    """
    image_pair, image_number, smooth_params, nuclei_params, edge_spot_params, perinuclear_params = (
        args
    )

    try:
        # Load images
        hoechst, miro = ImageLoader.load_image_pair(image_pair)

        # Apply smoothing (CellProfiler Smooth module) before segmentation
        if smooth_params is not None:
            hoechst = smooth_image(hoechst, smooth_params)

        # Segment nuclei - get both all-nuclei and interior-only versions
        # This fixes the border nuclei masking bug: we use ALL nuclei for masking
        # (so spots near border nuclei are properly excluded), but only interior
        # nuclei for Gini analysis (where truncated perinuclear rings are problematic)
        nuclei_labels_all, nuclei_labels_interior = segment_nuclei_keep_border(
            hoechst, nuclei_params
        )
        n_nuclei_all = nuclei_labels_all.max()
        n_nuclei_interior = nuclei_labels_interior.max()

        if n_nuclei_all == 0:
            metadata = image_pair.to_metadata_dict(image_number)
            return {
                "image_id": image_pair.image_id,
                "measurements": {
                    "Nuclei": pd.DataFrame(),
                    "Expand_Nuclei": pd.DataFrame(),
                    "Perinuclear_region": pd.DataFrame(),
                    "edge_spots": pd.DataFrame(),
                    "Image": pd.DataFrame([metadata]),
                },
            }

        # Create perinuclear regions using ALL nuclei (including border)
        # This ensures spots near border nuclei are masked out
        expand_nuclei_10px_all, expand_nuclei_15px_all, _ = create_perinuclear_regions(
            nuclei_labels_all, perinuclear_params
        )

        # Mask peripheral regions using ALL nuclei expansion
        # Centrosomal spots near border nuclei are now properly masked
        masked_miro = mask_peripheral_regions(miro, expand_nuclei_15px_all)

        # Detect edge spots on the properly masked MIRO
        edge_spots_labels = detect_edge_spots(masked_miro, edge_spot_params)
        n_edge_spots = edge_spots_labels.max()

        # Filter edge spots (Module 13)
        # CellProfiler measures and filters edge_spots on masked MIRO (Subtract_perinucleus_Miro160mer)
        if n_edge_spots > 0:
            edge_spots_labels = filter_edge_spots_by_edge_intensity(
                edge_spots_labels, masked_miro, min_edge_intensity=0.0, max_edge_intensity=1.0
            )

            # Filter by nuclei proximity - remove spots too close to perinuclear region
            edge_spots_labels = filter_edge_spots_by_nuclei_proximity(
                edge_spots_labels,
                expand_nuclei_15px_all,
                expansion_radius=3,
            )

        # Create perinuclear regions for INTERIOR nuclei only (for Gini analysis)
        # Truncated perinuclear rings at borders are problematic for Gini calculation
        expand_nuclei_10px, _, perinuclear_ring = create_perinuclear_regions(
            nuclei_labels_interior, perinuclear_params
        )

        # Measure properties using interior nuclei for Gini/dispersion metrics
        metadata = image_pair.to_metadata_dict(image_number)
        measurements = combine_measurements_for_export(
            nuclei_labels=nuclei_labels_interior,
            expand_nuclei_labels=expand_nuclei_10px,
            perinuclear_region_labels=perinuclear_ring,
            edge_spots_labels=edge_spots_labels,
            hoechst_image=hoechst,
            miro_image=miro,
            masked_miro_image=masked_miro,  # For edge_spots measurements (Module 12)
            metadata=metadata,
            nuclei_count_all=n_nuclei_all,
            nuclei_count_interior=n_nuclei_interior,
        )

        return {
            "image_id": image_pair.image_id,
            "measurements": measurements,
            "n_nuclei": n_nuclei_all,
            "n_edge_spots": edge_spots_labels.max(),
        }

    except Exception as e:
        return {"image_id": image_pair.image_id, "error": str(e)}


class Pipeline:
    """Complete image analysis pipeline."""

    def __init__(
        self,
        smooth_params: SmoothParams | None = None,
        nuclei_params: NucleiSegmentationParams | None = None,
        edge_spot_params: EdgeSpotParams | None = None,
        perinuclear_params: PerinuclearRegionParams | None = None,
    ):
        """
        Initialize pipeline with parameters.

        Args:
            smooth_params: Image smoothing parameters (None to disable smoothing)
            nuclei_params: Nuclei segmentation parameters (uses defaults if None)
            edge_spot_params: Edge spot detection parameters (uses defaults if None)
            perinuclear_params: Perinuclear region parameters (uses defaults if None)
        """
        self.smooth_params = smooth_params if smooth_params is not None else SmoothParams()
        self.nuclei_params = nuclei_params or NucleiSegmentationParams()
        self.edge_spot_params = edge_spot_params or EdgeSpotParams()
        self.perinuclear_params = perinuclear_params or PerinuclearRegionParams()
        self._current_exclusions: set[tuple[str, int]] = set()  # (well, xy) pairs to exclude
        self._figure_inclusions: set[tuple[str, str]] | None = None  # (date, well) pairs to include

    def process_image_pair(
        self, image_pair: ImagePair, image_number: int
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

        # Apply smoothing (CellProfiler Smooth module) before segmentation
        if self.smooth_params is not None:
            hoechst = smooth_image(hoechst, self.smooth_params)

        # Step 1: Segment nuclei (Module 7)
        # Get both all-nuclei and interior-only versions to fix border masking bug
        logger.debug("  Segmenting nuclei...")
        nuclei_labels_all, nuclei_labels_interior = segment_nuclei_keep_border(
            hoechst, self.nuclei_params
        )
        n_nuclei_all = nuclei_labels_all.max()
        n_nuclei_interior = nuclei_labels_interior.max()
        logger.debug(f"  Found {n_nuclei_all} nuclei ({n_nuclei_interior} interior)")

        if n_nuclei_all == 0:
            logger.warning(f"  No nuclei found in {image_pair.image_id}, skipping")
            return self._empty_measurements(image_pair, image_number)

        # Step 2: Create perinuclear regions using ALL nuclei (including border)
        # This ensures spots near border nuclei are properly masked out
        logger.debug("  Creating perinuclear regions (all nuclei for masking)...")
        expand_nuclei_10px_all, expand_nuclei_15px_all, _ = create_perinuclear_regions(
            nuclei_labels_all, self.perinuclear_params
        )

        # Step 3: Mask peripheral regions using ALL nuclei (Module 10)
        # Centrosomal spots near border nuclei are now properly masked
        logger.debug("  Masking peripheral regions...")
        masked_miro = mask_peripheral_regions(miro, expand_nuclei_15px_all)

        # Step 4: Detect edge spots (Module 11)
        logger.debug("  Detecting edge spots...")
        edge_spots_labels = detect_edge_spots(masked_miro, self.edge_spot_params)
        n_edge_spots = edge_spots_labels.max()
        logger.debug(f"  Found {n_edge_spots} edge spots")

        # Step 5: Filter edge spots by edge intensity (Module 13)
        # CellProfiler measures and filters edge_spots on masked MIRO (Subtract_perinucleus_Miro160mer)
        if n_edge_spots > 0:
            logger.debug("  Filtering edge spots by edge intensity...")
            edge_spots_labels = filter_edge_spots_by_edge_intensity(
                edge_spots_labels,
                masked_miro,  # Use masked image per CellProfiler Module 12-13
                min_edge_intensity=0.0,
                max_edge_intensity=1.0,
            )
            n_edge_spots_filtered = edge_spots_labels.max()
            logger.debug(f"  {n_edge_spots_filtered} spots after edge intensity filter")

            # Filter by nuclei proximity - remove spots too close to perinuclear region
            logger.debug("  Filtering edge spots by nuclei proximity...")
            edge_spots_labels = filter_edge_spots_by_nuclei_proximity(
                edge_spots_labels,
                expand_nuclei_15px_all,
                expansion_radius=3,
            )
            n_edge_spots_proximity = edge_spots_labels.max()
            logger.debug(f"  {n_edge_spots_proximity} spots after nuclei proximity filter")

        # Step 6: Create perinuclear regions for INTERIOR nuclei only (for Gini)
        # Truncated perinuclear rings at borders are problematic for Gini calculation
        logger.debug("  Creating perinuclear regions (interior nuclei for Gini)...")
        expand_nuclei_10px, _, perinuclear_ring = create_perinuclear_regions(
            nuclei_labels_interior, self.perinuclear_params
        )

        # Step 7: Measure all properties (Modules 12-20)
        logger.debug("  Measuring properties...")
        metadata = image_pair.to_metadata_dict(image_number)

        measurements = combine_measurements_for_export(
            nuclei_labels=nuclei_labels_interior,
            expand_nuclei_labels=expand_nuclei_10px,
            perinuclear_region_labels=perinuclear_ring,
            edge_spots_labels=edge_spots_labels,
            hoechst_image=hoechst,
            miro_image=miro,
            masked_miro_image=masked_miro,  # For edge_spots measurements (Module 12)
            metadata=metadata,
            nuclei_count_all=n_nuclei_all,
            nuclei_count_interior=n_nuclei_interior,
        )

        logger.info(f"  Complete: {n_nuclei_all} nuclei, {edge_spots_labels.max()} edge spots")

        return measurements

    def _empty_measurements(
        self, image_pair: ImagePair, image_number: int
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
            "Nuclei": pd.DataFrame(),
            "Expand_Nuclei": pd.DataFrame(),
            "Perinuclear_region": pd.DataFrame(),
            "edge_spots": pd.DataFrame(),
            "Image": pd.DataFrame([metadata]),
        }

    def process_date_batch(
        self, input_dir: Path, output_dir: Path, date: str, n_workers: int = 1
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

        # Filter by figure inclusions (--figures-only mode)
        if self._figure_inclusions is not None:
            original_count = len(image_pairs)
            image_pairs = [
                ip for ip in image_pairs if should_include_image(ip, self._figure_inclusions)
            ]
            filtered_count = original_count - len(image_pairs)
            if filtered_count > 0:
                logger.info(
                    f"Filtered {filtered_count} images not in figure sheets "
                    f"({len(image_pairs)} remaining)"
                )

        # Filter out excluded images
        if self._current_exclusions:
            original_count = len(image_pairs)
            image_pairs = [
                ip for ip in image_pairs if not should_exclude_image(ip, self._current_exclusions)
            ]
            excluded_count = original_count - len(image_pairs)
            if excluded_count > 0:
                logger.info(
                    f"Excluded {excluded_count} images based on config: {self._current_exclusions}"
                )

        # Prepare args for parallel processing
        args_list = [
            (
                pair,
                i,
                self.smooth_params,
                self.nuclei_params,
                self.edge_spot_params,
                self.perinuclear_params,
            )
            for i, pair in enumerate(image_pairs, start=1)
        ]

        all_measurements = []

        if n_workers == 1:
            # Sequential processing
            for args in args_list:
                result = _process_single_image(args)
                if "error" in result:
                    logger.error(f"Error processing {result['image_id']}: {result['error']}")
                else:
                    logger.info(
                        f"  {result['image_id']}: {result.get('n_nuclei', 0)} nuclei, {result.get('n_edge_spots', 0)} edge spots"
                    )
                    all_measurements.append(result["measurements"])
        else:
            # Parallel processing
            logger.info(f"Using {n_workers} parallel workers")
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(_process_single_image, args): args[0].image_id
                    for args in args_list
                }

                for future in as_completed(futures):
                    image_id = futures[future]
                    try:
                        result = future.result()
                        if "error" in result:
                            logger.error(
                                f"Error processing {result['image_id']}: {result['error']}"
                            )
                        else:
                            logger.info(
                                f"  {result['image_id']}: {result.get('n_nuclei', 0)} nuclei, {result.get('n_edge_spots', 0)} edge spots"
                            )
                            all_measurements.append(result["measurements"])
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
        config: dict[str, dict[str, Any]] | None = None,
        n_workers: int = 1,
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
        if "otsu_correction_factor" in params:
            self.nuclei_params.otsu_correction_factor = params["otsu_correction_factor"]
        if "diameter_min" in params:
            self.nuclei_params.diameter_min = params["diameter_min"]
        if "diameter_max" in params:
            self.nuclei_params.diameter_max = params["diameter_max"]
        if "min_distance" in params:
            self.nuclei_params.min_distance = params["min_distance"]
        if "threshold_smoothing_scale" in params:
            self.nuclei_params.threshold_smoothing_scale = params["threshold_smoothing_scale"]
        if "threshold_lower_bound" in params:
            self.nuclei_params.threshold_lower_bound = params["threshold_lower_bound"]
        if "threshold_upper_bound" in params:
            self.nuclei_params.threshold_upper_bound = params["threshold_upper_bound"]

        # Update smooth params (CellProfiler Smooth module)
        if "smooth_artifact_diameter" in params:
            if self.smooth_params is None:
                self.smooth_params = SmoothParams()
            self.smooth_params.artifact_diameter = params["smooth_artifact_diameter"]
        if "smooth_method" in params:
            if self.smooth_params is None:
                self.smooth_params = SmoothParams()
            self.smooth_params.method = params["smooth_method"]

        # Update edge spot params
        if "edge_spot_diameter_min" in params:
            self.edge_spot_params.diameter_min = params["edge_spot_diameter_min"]
        if "edge_spot_diameter_max" in params:
            self.edge_spot_params.diameter_max = params["edge_spot_diameter_max"]
        if "edge_spot_correction_factor" in params:
            self.edge_spot_params.correction_factor = params["edge_spot_correction_factor"]

        # Update perinuclear region params
        if "inner_expansion" in params:
            self.perinuclear_params.inner_expansion = params["inner_expansion"]
        if "outer_expansion" in params:
            self.perinuclear_params.outer_expansion = params["outer_expansion"]

        # Update exclusions for current date
        self._current_exclusions = params.get("exclusions", set())

    def load_figure_inclusions(self, config_path: Path) -> None:
        """
        Load figure inclusions from config file.

        When set, only images whose (date, well) appear in Figure sheets
        will be processed.

        Args:
            config_path: Path to Excel config file with Figure sheets
        """
        self._figure_inclusions = extract_figure_inclusions(config_path)
        logger.info(f"Loaded {len(self._figure_inclusions)} figure inclusions")


def main():
    """Main entry point for CLI."""
    # Configure logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    parser = argparse.ArgumentParser(
        description="Run CellProfiler pipeline in Python",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all dates with default parameters
  edge-spot-pipeline --input inputs/ --output results/

  # Process specific date
  edge-spot-pipeline --input inputs/ --output results/ --date 20231115

  # Use custom parameter configuration
  edge-spot-pipeline --input inputs/ --output results/ --config params.json

  # Skip aggregation (only run image processing)
  edge-spot-pipeline --input inputs/ --output results/ --skip-aggregate

  # Enable time-course mode for aggregation
  edge-spot-pipeline --input inputs/ --output results/ --t-varies

  # Verbose logging
  edge-spot-pipeline --input inputs/ --output results/ --verbose
        """,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input directory containing date subdirectories with TIF files",
    )

    parser.add_argument(
        "--output", "-o", type=Path, required=True, help="Output directory for CSV files"
    )

    parser.add_argument("--date", "-d", type=str, help="Process only this specific date (optional)")

    parser.add_argument(
        "--config", "-c", type=Path, help="Path to parameter configuration file (JSON or YAML)"
    )

    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    parser.add_argument(
        "--workers",
        "-w",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1, use -1 for all CPUs)",
    )

    parser.add_argument(
        "--skip-aggregate",
        action="store_true",
        help="Skip aggregation step (only run image processing)",
    )

    parser.add_argument(
        "--t-varies", action="store_true", help="Enable time-course mode for aggregation"
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots during aggregation (requires matplotlib)",
    )

    parser.add_argument(
        "--figures-only",
        action="store_true",
        help="Only process images that appear in Figure sheets (requires --config)",
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

    # Load figure inclusions if --figures-only is set
    if args.figures_only:
        if not args.config:
            logger.error("--figures-only requires --config to specify the Excel config file")
            return
        pipeline.load_figure_inclusions(args.config)

    # Process images
    start_time = time.time()

    if args.date:
        # Process single date - apply config if provided
        if config:
            date_params = get_parameters_for_date(config, args.date)
            pipeline._update_parameters(date_params)
            logger.info(f"Using parameters for date {args.date}: {date_params}")
        pipeline.process_date_batch(args.input, args.output, args.date, n_workers=n_workers)
    else:
        # Process all dates
        pipeline.process_all_dates(args.input, args.output, config, n_workers=n_workers)

    elapsed_time = time.time() - start_time
    logger.info(f"Pipeline complete in {elapsed_time:.2f} seconds")

    # Run aggregation unless skipped
    if not args.skip_aggregate:
        from edge_spot_analyser.aggregation import AggregationConfig, Aggregator

        logger.info("Running aggregation...")
        agg_config = AggregationConfig(t_varies=args.t_varies, plot=args.plot)
        aggregator = Aggregator(agg_config)

        output_path = Path(args.output)
        if args.date:
            # Single date
            date_output = output_path / args.date
            if date_output.exists():
                aggregator.process_date(date_output)
        else:
            # All dates
            for date_dir in output_path.iterdir():
                if date_dir.is_dir():
                    aggregator.process_date(date_dir)

        logger.info("Aggregation complete")


if __name__ == "__main__":
    main()
