"""
Edge Spot Analyser - Python port of CellProfiler pipeline.

This package analyzes mitochondrial distribution around nuclei in microscopy images.
"""

__version__ = "1.0.0"

from edge_spot_analyser.aggregation import (
    AggregationConfig,
    Aggregator,
)
from edge_spot_analyser.figure_aggregation import FigureAggregator
from edge_spot_analyser.io_utils import (
    CSVExporter,
    FileDiscovery,
    ImageLoader,
    ImagePair,
    get_parameters_for_date,
    load_parameter_config,
    parse_exclusions,
    should_exclude_image,
)
from edge_spot_analyser.measurements import (
    ImageMeasurements,
    IntensityMeasurements,
    combine_measurements_for_export,
    measure_all_object_properties,
)
from edge_spot_analyser.pipeline import Pipeline
from edge_spot_analyser.segmentation import (
    EdgeSpotParams,
    NucleiSegmentationParams,
    PerinuclearRegionParams,
    SmoothParams,
    create_perinuclear_regions,
    detect_edge_spots,
    filter_edge_spots_by_edge_intensity,
    mask_peripheral_regions,
    segment_nuclei,
    smooth_image,
)

__all__ = [
    "__version__",
    # Segmentation
    "SmoothParams",
    "smooth_image",
    "NucleiSegmentationParams",
    "EdgeSpotParams",
    "PerinuclearRegionParams",
    "segment_nuclei",
    "create_perinuclear_regions",
    "mask_peripheral_regions",
    "detect_edge_spots",
    "filter_edge_spots_by_edge_intensity",
    # Measurements
    "IntensityMeasurements",
    "ImageMeasurements",
    "measure_all_object_properties",
    "combine_measurements_for_export",
    # I/O
    "ImagePair",
    "FileDiscovery",
    "ImageLoader",
    "CSVExporter",
    "load_parameter_config",
    "get_parameters_for_date",
    "parse_exclusions",
    "should_exclude_image",
    # Pipeline
    "Pipeline",
    # Aggregation
    "Aggregator",
    "AggregationConfig",
    "FigureAggregator",
]
