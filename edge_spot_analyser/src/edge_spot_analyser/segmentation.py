"""
Segmentation module for the CellProfiler pipeline port.

This module implements the core segmentation steps from the CellProfiler pipeline:
- Smooth module: Image smoothing/noise reduction (median filter)
- Module 7: Nuclei segmentation using 3-class Otsu thresholding
- Module 8-10, 14: Perinuclear region creation
- Module 11: Edge spot detection using Robust Background thresholding
- Module 13: Edge spot filtering by edge intensity

Author: Generated from CellProfiler pipeline 231115_combined_pipeline_new_nomito_fixed_moments
Pipeline Version: v6
"""

import logging
from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter, median_filter
from skimage import filters, measure, morphology, segmentation
from skimage.feature import peak_local_max

logger = logging.getLogger(__name__)


@dataclass
class SmoothParams:
    """Parameters for image smoothing (CellProfiler Smooth module)."""

    method: str = "median"  # Smoothing method: "median" or "gaussian"
    artifact_diameter: int = 3  # Size of artifacts to remove (for median filter)
    sigma: float = 1.0  # Sigma for Gaussian smoothing (if method="gaussian")


def smooth_image(image: np.ndarray, params: SmoothParams | None = None) -> np.ndarray:
    """
    Apply smoothing to an image for noise reduction.

    Implements CellProfiler's Smooth module with median filter option.
    This is a preprocessing step applied BEFORE segmentation.

    Args:
        image: 2D float array (0-1 normalized)
        params: Smoothing parameters (uses defaults if None)

    Returns:
        Smoothed image

    References:
        CellProfiler Smooth module with "Median Filter" method
    """
    if params is None:
        params = SmoothParams()

    if params.method == "median":
        if params.artifact_diameter <= 0:
            return image
        # CellProfiler uses a disk-shaped structuring element
        # artifact_diameter is the diameter, so radius = diameter // 2
        radius = max(1, params.artifact_diameter // 2)
        footprint = morphology.disk(radius)
        return median_filter(image, footprint=footprint)
    elif params.method == "gaussian":
        return gaussian_filter(image, sigma=params.sigma)
    else:
        raise ValueError(f"Unknown smoothing method: {params.method}")


@dataclass
class NucleiSegmentationParams:
    """Parameters for nuclei segmentation (Module 7)."""

    diameter_min: int = 30  # Minimum nucleus diameter in pixels
    diameter_max: int = 100  # Maximum nucleus diameter in pixels
    otsu_correction_factor: float = 0.45  # Correction factor for Otsu threshold
    threshold_lower_bound: float = 0.00195  # Absolute minimum threshold
    threshold_upper_bound: float = 0.01  # Absolute maximum threshold
    min_distance: int = 30  # Minimum distance between nuclei for declumping
    threshold_smoothing_scale: float = (
        2.0  # CP "Threshold smoothing scale" - Gaussian smoothing before threshold
    )
    declump_smoothing_filter: int = 10  # CP "Size of smoothing filter" - for declumping
    discard_border_objects: bool = True  # Remove objects touching image border

    @property
    def size_min(self) -> int:
        """Minimum nucleus area in pixels."""
        return int(np.pi * (self.diameter_min / 2) ** 2)

    @property
    def size_max(self) -> int:
        """Maximum nucleus area in pixels."""
        return int(np.pi * (self.diameter_max / 2) ** 2)


@dataclass
class EdgeSpotParams:
    """Parameters for edge spot detection (Module 11)."""

    diameter_min: int = 5  # Minimum spot diameter in pixels
    diameter_max: int = 80  # Maximum spot diameter in pixels
    lower_outlier_fraction: float = 0.3  # Fraction of lowest pixels to remove
    upper_outlier_fraction: float = 0.1  # Fraction of highest pixels to remove
    num_deviations: float = 6.0  # Number of standard deviations above robust mean
    correction_factor: float = 2.0  # Multiplier applied to calculated threshold
    threshold_lower_bound: float = 0.0035  # Absolute minimum threshold
    threshold_smoothing: float = 1.3  # Gaussian smoothing applied to image before thresholding
    discard_border_objects: bool = True  # Remove objects touching image border

    @property
    def size_min(self) -> int:
        """Minimum spot area in pixels."""
        return int(np.pi * (self.diameter_min / 2) ** 2)

    @property
    def size_max(self) -> int:
        """Maximum spot area in pixels."""
        return int(np.pi * (self.diameter_max / 2) ** 2)


@dataclass
class PerinuclearRegionParams:
    """Parameters for perinuclear region creation (Modules 8-10, 14)."""

    inner_expansion: int = 10  # Pixels to expand for inner boundary
    outer_expansion: int = 15  # Pixels to expand for outer boundary (masking)


def segment_nuclei(
    hoechst_image: np.ndarray, params: NucleiSegmentationParams | None = None
) -> np.ndarray:
    """
    Segment nuclei from Hoechst channel using 3-class Otsu thresholding.

    Implements CellProfiler Module 7 (IdentifyPrimaryObjects → Nuclei).

    Note: For noise reduction, apply smooth_image() to the Hoechst image
    BEFORE calling this function (matching CellProfiler's Smooth module).

    Algorithm:
    1. Apply Gaussian smoothing to image (σ = threshold_smoothing_scale/2)
    2. Calculate 3-class Otsu threshold (uses upper threshold)
    3. Apply correction factor (default 0.45) and clip to bounds
    4. Create binary mask
    5. Declump using intensity-based watershed
    6. Filter by size (30-100 pixel diameter)
    7. Remove border objects

    Args:
        hoechst_image: 2D float array (0-1 normalized) of Hoechst/nuclear channel
        params: Segmentation parameters (uses defaults if None)

    Returns:
        Labeled image where each nucleus has a unique integer label (background=0)

    References:
        CellProfiler IdentifyPrimaryObjects with 3-class Otsu thresholding
    """
    if params is None:
        params = NucleiSegmentationParams()

    # Step 1: Apply Gaussian smoothing to image before thresholding
    # Note: Median filter preprocessing (CellProfiler Smooth module) should be
    # applied BEFORE calling this function using smooth_image()
    smoothed = hoechst_image
    # CP "Threshold smoothing scale" of 2.0 → sigma = 1.0 (scale/2)
    sigma = params.threshold_smoothing_scale / 2.0
    smoothed = gaussian_filter(smoothed, sigma=sigma)

    # Step 2: Clip bright artifacts before thresholding
    # Some images have very bright dead cells/artifacts that skew the Otsu threshold.
    # If max intensity >> 99th percentile, clip to remove these outliers.
    p99 = np.percentile(smoothed, 99)
    if smoothed.max() > p99 * 3:
        logger.debug(f"Clipping bright artifacts: max={smoothed.max():.4f}, p99={p99:.4f}")
        smoothed = np.clip(smoothed, 0, p99 * 1.5)

    # Step 3: Calculate 3-class Otsu thresholds using skimage
    # threshold_multiotsu returns array of thresholds dividing into N classes
    # For 3 classes with middle→background, we use the upper threshold (index 1)
    # CellProfiler uses 128 histogram bins by default
    try:
        thresholds = filters.threshold_multiotsu(smoothed, classes=3, nbins=128)
        threshold = thresholds[1]  # Upper threshold (middle class → background)
    except (ValueError, IndexError):
        # Fallback if image has insufficient dynamic range
        logger.warning("3-class Otsu failed (low dynamic range), falling back to 2-class Otsu")
        threshold = filters.threshold_otsu(smoothed)

    # Step 4: Apply correction factor and clip to bounds
    threshold *= params.otsu_correction_factor
    threshold = np.clip(threshold, params.threshold_lower_bound, params.threshold_upper_bound)

    # Step 5: Create binary mask
    binary = smoothed > threshold

    # Step 6: Declump using intensity-based watershed
    # Note: Fill holes AFTER declumping (CellProfiler: "After declumping only")
    # CellProfiler auto smoothing: sigma = min_diameter / 2.35 / 2
    declump_smoothing_sigma = params.diameter_min / 4.7
    labels = _declump_objects_watershed(
        binary,
        intensity_image=hoechst_image,
        min_distance=params.min_distance,
        smoothing_sigma=declump_smoothing_sigma,
    )

    # Step 7: Fill holes in each labeled object (after declumping)
    labels = _fill_holes_in_labels(labels)

    # Step 8: Filter by size
    labels = _filter_objects_by_size(labels, params.size_min, params.size_max)

    # Step 9: Remove border objects
    if params.discard_border_objects:
        labels = segmentation.clear_border(labels)

    # Relabel consecutively
    labels = measure.label(labels > 0)

    return labels


def _declump_objects_watershed(
    binary_mask: np.ndarray,
    intensity_image: np.ndarray,
    min_distance: int,
    smoothing_sigma: float | None = None,
) -> np.ndarray:
    """
    Declump touching objects using intensity-based watershed.

    Args:
        binary_mask: Binary mask of objects
        intensity_image: Original intensity image for watershed
        min_distance: Minimum distance between object centers
        smoothing_sigma: Gaussian smoothing sigma for intensity image before watershed.
            If None, no smoothing is applied.

    Returns:
        Labeled image with separated objects
    """
    # Compute distance transform
    distance = ndimage.distance_transform_edt(binary_mask)

    # Find local maxima as seeds (use min_distance to prevent over-segmentation)
    coords = morphology.local_maxima(distance, indices=True)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords)] = True

    # Apply minimum distance constraint
    # Dilate each maximum by min_distance and keep only non-overlapping ones
    coords_array = peak_local_max(distance, min_distance=min_distance, labels=binary_mask)

    # Create markers from peaks
    markers = np.zeros(binary_mask.shape, dtype=int)
    for i, coord in enumerate(coords_array, start=1):
        markers[tuple(coord)] = i

    # Smooth intensity image before watershed (CellProfiler "Size of smoothing filter")
    if smoothing_sigma is not None and smoothing_sigma > 0:
        intensity_for_watershed = gaussian_filter(intensity_image, sigma=smoothing_sigma)
    else:
        intensity_for_watershed = intensity_image

    # Watershed using inverted intensity (watershed finds minima)
    # CellProfiler uses "intensity-based" watershed, meaning watershed on inverted intensity
    labels = segmentation.watershed(-intensity_for_watershed, markers, mask=binary_mask)

    return labels


def _fill_holes_in_labels(labels: np.ndarray) -> np.ndarray:
    """
    Fill holes in each labeled object independently.

    This matches CellProfiler's "Fill holes after declumping only" behavior.

    Args:
        labels: Labeled image where each object has unique integer label

    Returns:
        Labeled image with holes filled in each object
    """
    output = labels.copy()  # Start with original labels to preserve existing assignments

    for label_id in range(1, labels.max() + 1):
        # Extract single object mask
        obj_mask = labels == label_id
        # Fill holes in this object
        filled = ndimage.binary_fill_holes(obj_mask)
        # Only fill holes that don't overlap with other objects
        # This prevents merging adjacent objects whose filled holes would touch
        new_pixels = filled & (output == 0)
        output[new_pixels] = label_id

    return output


def _filter_objects_by_size(labels: np.ndarray, size_min: int, size_max: int) -> np.ndarray:
    """
    Filter labeled objects by area.

    Args:
        labels: Labeled image
        size_min: Minimum object area in pixels
        size_max: Maximum object area in pixels

    Returns:
        Filtered labeled image
    """
    # Get region properties
    props = measure.regionprops(labels)

    # Create mask of objects to keep
    keep_labels = [prop.label for prop in props if size_min <= prop.area <= size_max]

    # Create filtered image
    filtered = np.zeros_like(labels)
    for label in keep_labels:
        filtered[labels == label] = label

    return filtered


def create_perinuclear_regions(
    nuclei_labels: np.ndarray, params: PerinuclearRegionParams | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create perinuclear regions by expanding nuclei.

    Implements CellProfiler Modules 8-10, 14:
    - Module 8: Expand nuclei by 10 pixels (Expand_Nuclei)
    - Module 9: Expand nuclei by 15 pixels (Expand_Nuclei_for_mask, for masking)
    - Module 14: Create tertiary object (Perinuclear_region)
      - Larger object: Expand_Nuclei (10px expansion)
      - Smaller object: Nuclei (original, eroded by 1px)
      - Result: Ring between eroded nuclei and 10px expansion

    Args:
        nuclei_labels: Labeled image of nuclei
        params: Perinuclear region parameters (uses defaults if None)

    Returns:
        Tuple of (expand_nuclei_10px, expand_nuclei_15px, perinuclear_ring)
        - expand_nuclei_10px: Nuclei expanded by 10 pixels
        - expand_nuclei_15px: Nuclei expanded by 15 pixels
        - perinuclear_ring: Ring between eroded nuclei and 10px expansion
    """
    if params is None:
        params = PerinuclearRegionParams()

    # Module 8: Expand by 10 pixels
    expand_nuclei_10px = morphology.dilation(nuclei_labels, morphology.disk(params.inner_expansion))

    # Module 9: Expand by 15 pixels (for masking peripheral regions)
    expand_nuclei_15px = morphology.dilation(nuclei_labels, morphology.disk(params.outer_expansion))

    # Module 14: Create tertiary object (perinuclear ring)
    # CellProfiler IdentifyTertiaryObjects with:
    #   - Larger: Expand_Nuclei (10px expansion)
    #   - Smaller: Nuclei (original)
    #   - "Shrink smaller object prior to subtraction": Yes (erode by 1px)
    # Result: Ring between (nuclei eroded by 1px) and (10px expansion)
    eroded_nuclei = morphology.erosion(nuclei_labels, morphology.disk(1))
    perinuclear_ring = np.where(
        (expand_nuclei_10px > 0) & (eroded_nuclei == 0),
        expand_nuclei_10px,  # Keep the nucleus label from expansion
        0,
    )

    return expand_nuclei_10px, expand_nuclei_15px, perinuclear_ring


def mask_peripheral_regions(miro_image: np.ndarray, expand_nuclei_15px: np.ndarray) -> np.ndarray:
    """
    Mask out perinuclear regions from MIRO160mer image.

    Implements CellProfiler Module 10 (MaskImage → Subtract_perinucleus_Miro160mer).
    Creates an inverted mask that keeps only peripheral areas (outside 15px expansion).

    Args:
        miro_image: 2D float array (0-1 normalized) of MIRO160mer channel
        expand_nuclei_15px: Labeled image of nuclei expanded by 15 pixels

    Returns:
        Masked MIRO160mer image with perinuclear regions set to 0
    """
    # Create inverted mask: keep only regions outside 15px expansion
    peripheral_mask = expand_nuclei_15px == 0

    # Apply mask
    masked_miro = miro_image.copy()
    masked_miro[~peripheral_mask] = 0

    return masked_miro


def robust_background_threshold(
    image: np.ndarray,
    lower_outlier_fraction: float = 0.3,
    upper_outlier_fraction: float = 0.1,
    num_deviations: float = 6.0,
    correction_factor: float = 2.0,
    lower_bound: float = 0.0035,
    upper_bound: float = 1.0,
) -> float:
    """
    Calculate threshold using CellProfiler's Robust Background algorithm.

    This is a CellProfiler-specific thresholding method that's critical for
    detecting very bright mitochondrial spots.

    Algorithm:
    1. Extract non-zero pixels and sort
    2. Remove outliers: lowest 30%, highest 10%
    3. Calculate robust statistics on remaining 60%:
       - robust_mean = MEDIAN (not mean!)
       - robust_std = standard deviation
    4. Calculate threshold:
       threshold = robust_mean + (num_deviations × robust_std)
       threshold *= correction_factor
       threshold = clip(threshold, lower_bound, upper_bound)

    Args:
        image: 2D float array (0-1 normalized)
        lower_outlier_fraction: Fraction of lowest pixels to remove (default 0.3)
        upper_outlier_fraction: Fraction of highest pixels to remove (default 0.1)
        num_deviations: Number of standard deviations above robust mean (default 6.0)
        correction_factor: Multiplier for threshold (default 2.0)
        lower_bound: Absolute minimum threshold (default 0.0035)
        upper_bound: Absolute maximum threshold (default 1.0)

    Returns:
        Calculated threshold value

    References:
        CellProfiler Robust Background thresholding method
        See MODULE_11_DEFINITIVE.md for detailed algorithm
    """
    # Step 1: Extract non-zero pixels and sort
    pixels = image[image > 0].flatten()
    if len(pixels) == 0:
        return lower_bound

    pixels_sorted = np.sort(pixels)
    n_pixels = len(pixels_sorted)

    # Step 2: Remove outliers
    lower_cutoff = int(n_pixels * lower_outlier_fraction)
    upper_cutoff = int(n_pixels * (1 - upper_outlier_fraction))

    robust_pixels = pixels_sorted[lower_cutoff:upper_cutoff]

    if len(robust_pixels) == 0:
        return lower_bound

    # Step 3: Calculate robust statistics
    # CRITICAL: Use MEDIAN as robust_mean, not mean!
    robust_mean = np.median(robust_pixels)
    robust_std = np.std(robust_pixels)

    # Step 4: Calculate threshold
    threshold = robust_mean + (num_deviations * robust_std)
    threshold *= correction_factor
    threshold = np.clip(threshold, lower_bound, upper_bound)

    return threshold


def detect_edge_spots(masked_miro: np.ndarray, params: EdgeSpotParams | None = None) -> np.ndarray:
    """
    Detect bright peripheral mitochondrial spots using Robust Background thresholding.

    Implements CellProfiler Module 11 (IdentifyPrimaryObjects → edge_spots).

    Algorithm:
    1. Apply Gaussian smoothing to image (σ = threshold_smoothing / 2.0)
    2. Calculate threshold using Robust Background method on smoothed image
    3. Create binary mask: smoothed > threshold
    4. Label connected components (NO watershed declumping)
    5. Fill holes
    6. Filter by size (5-80 pixel diameter)
    7. Remove border objects

    Args:
        masked_miro: 2D float array of masked MIRO160mer (peripheral regions only)
        params: Edge spot detection parameters (uses defaults if None)

    Returns:
        Labeled image where each spot has a unique integer label (background=0)

    References:
        CellProfiler IdentifyPrimaryObjects with Robust Background thresholding
        See MODULE_11_DEFINITIVE.md for complete algorithm
    """
    if params is None:
        params = EdgeSpotParams()

    # Step 1: Apply Gaussian smoothing to image before thresholding (CellProfiler behavior)
    # CellProfiler "Threshold smoothing scale" of N → sigma = N/2 (matches nuclei segmentation)
    sigma = params.threshold_smoothing / 2.0
    smoothed = gaussian_filter(masked_miro, sigma=sigma)

    # Step 2: Calculate threshold using Robust Background on smoothed image
    threshold = robust_background_threshold(
        smoothed,
        lower_outlier_fraction=params.lower_outlier_fraction,
        upper_outlier_fraction=params.upper_outlier_fraction,
        num_deviations=params.num_deviations,
        correction_factor=params.correction_factor,
        lower_bound=params.threshold_lower_bound,
        upper_bound=1.0,
    )

    # Step 3: Create binary mask from smoothed image
    binary = smoothed > threshold

    # Step 4: Label connected components
    # IMPORTANT: Module 11 uses "Method to distinguish clumped objects: None"
    # This means NO watershed, just simple connected component labeling
    labels = measure.label(binary)

    # Step 6: Fill holes in each object
    labels = _fill_holes_in_labels(labels)

    # Step 7: Filter by size
    labels = _filter_objects_by_size(labels, params.size_min, params.size_max)

    # Step 8: Remove border objects
    if params.discard_border_objects:
        labels = segmentation.clear_border(labels)

    # Relabel consecutively
    labels = measure.label(labels > 0)

    return labels


def filter_edge_spots_by_edge_intensity(
    edge_spots_labels: np.ndarray,
    miro_image: np.ndarray,
    min_edge_intensity: float = 0.0,
    max_edge_intensity: float = 1.0,
) -> np.ndarray:
    """
    Filter edge spots by their edge intensity.

    Implements CellProfiler Module 13 (FilterObjects).

    Edge intensity is calculated by:
    1. Dilating each spot by 1 pixel
    2. Measuring mean intensity in the border region (dilated - original)

    Args:
        edge_spots_labels: Labeled image of edge spots
        miro_image: Original MIRO160mer image
        min_edge_intensity: Minimum edge intensity (default 0.0)
        max_edge_intensity: Maximum edge intensity (default 1.0)

    Returns:
        Filtered labeled image
    """
    if edge_spots_labels.max() == 0:
        return edge_spots_labels

    # Calculate edge intensity for each spot
    filtered = np.zeros_like(edge_spots_labels)

    for region in measure.regionprops(edge_spots_labels):
        # Get object mask
        mask = edge_spots_labels == region.label

        # Dilate by 1 pixel
        dilated_mask = morphology.binary_dilation(mask, morphology.disk(1))

        # Edge is dilated - original
        edge_mask = dilated_mask & ~mask

        # Calculate mean edge intensity
        if edge_mask.sum() > 0:
            edge_intensity = miro_image[edge_mask].mean()
        else:
            edge_intensity = 0.0

        # Keep if within range
        if min_edge_intensity <= edge_intensity <= max_edge_intensity:
            filtered[mask] = region.label

    # Relabel consecutively
    filtered = measure.label(filtered > 0)

    return filtered
