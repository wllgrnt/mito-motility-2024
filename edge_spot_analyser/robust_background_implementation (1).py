"""
Analysis of CellProfiler's Robust Background thresholding method

Based on CellProfiler source code:
https://github.com/CellProfiler/CellProfiler/blob/master/src/subpackages/core/cellprofiler_core/threshold/_robust_background_threshold.py
"""

import numpy as np
from scipy import ndimage
from skimage import filters


def robust_background_threshold(image,
                                 lower_outlier_fraction=0.3,
                                 upper_outlier_fraction=0.1,
                                 averaging_method='median',
                                 variance_method='std',
                                 num_deviations=6,
                                 correction_factor=1.0):
    """
    CellProfiler's Robust Background thresholding algorithm

    This method calculates a threshold by:
    1. Removing outliers from the image histogram
    2. Computing robust statistics (mean/median and std)
    3. Setting threshold = robust_mean + (num_deviations * robust_std)
    4. Applying correction factor

    Args:
        image: Input image (grayscale, 0-1 range)
        lower_outlier_fraction: Fraction of low intensities to exclude (0.0-1.0)
        upper_outlier_fraction: Fraction of high intensities to exclude (0.0-1.0)
        averaging_method: 'mean' or 'median'
        variance_method: 'std' or 'mad' (median absolute deviation)
        num_deviations: Number of standard deviations above mean
        correction_factor: Multiplier for final threshold

    Returns:
        threshold: Single threshold value
    """
    # Flatten and sort image
    pixels = image.flatten()
    pixels = pixels[pixels > 0]  # Ignore zero pixels (background)

    if len(pixels) == 0:
        return 0.0

    pixels_sorted = np.sort(pixels)
    n_pixels = len(pixels_sorted)

    # Calculate outlier indices
    lower_idx = int(n_pixels * lower_outlier_fraction)
    upper_idx = int(n_pixels * (1.0 - upper_outlier_fraction))

    # Get robust pixels (exclude outliers)
    robust_pixels = pixels_sorted[lower_idx:upper_idx]

    if len(robust_pixels) == 0:
        return 0.0

    # Calculate robust mean
    if averaging_method.lower() == 'median':
        robust_mean = np.median(robust_pixels)
    else:  # mean
        robust_mean = np.mean(robust_pixels)

    # Calculate robust variance
    if variance_method.lower() == 'mad':
        # Median Absolute Deviation
        deviations = np.abs(robust_pixels - np.median(robust_pixels))
        robust_std = np.median(deviations) * 1.4826  # Scale to match std
    else:  # std
        robust_std = np.std(robust_pixels, ddof=1)

    # Calculate threshold
    threshold = robust_mean + (num_deviations * robust_std)

    # Apply correction factor
    threshold *= correction_factor

    return threshold


def identify_edge_spots_correct(masked_image,
                                min_diameter=5,
                                max_diameter=80,
                                threshold_smoothing_scale=1.3,
                                threshold_correction_factor=2.0,
                                lower_outlier_fraction=0.3,
                                upper_outlier_fraction=0.1,
                                num_deviations=6,
                                lower_bound=0.0035,
                                upper_bound=1.0):
    """
    Complete implementation of edge_spots identification (Module 11)

    CRITICAL: This implementation matches CellProfiler's exact behavior
    """
    from skimage import measure, morphology, segmentation

    img = masked_image.copy()

    # Step 1: Calculate Robust Background threshold
    threshold = robust_background_threshold(
        img,
        lower_outlier_fraction=lower_outlier_fraction,
        upper_outlier_fraction=upper_outlier_fraction,
        averaging_method='median',
        variance_method='std',
        num_deviations=num_deviations,
        correction_factor=threshold_correction_factor
    )

    # Step 2: Apply bounds
    threshold = np.clip(threshold, lower_bound, upper_bound)

    print(f"Calculated threshold: {threshold:.6f}")

    # Step 3: Apply threshold
    binary = img > threshold

    # Step 4: Smooth the BINARY image (not the grayscale)
    # This is the "Threshold smoothing scale" parameter
    if threshold_smoothing_scale > 0:
        # Smooth the binary mask
        binary_float = binary.astype(float)
        binary_smoothed = filters.gaussian(binary_float, sigma=threshold_smoothing_scale)
        binary = binary_smoothed > 0.5

    # Step 5: Size filtering (before labeling)
    min_size = int(np.pi * (min_diameter/2)**2)
    max_size = int(np.pi * (max_diameter/2)**2)

    binary = morphology.remove_small_objects(binary, min_size=min_size)

    # Step 6: Label objects
    # Method to distinguish clumped objects: None
    # This means NO watershed, NO declumping - just label connected components
    labels = morphology.label(binary)

    # Step 7: Filter by size (remove large objects)
    props = measure.regionprops(labels)
    for prop in props:
        if prop.area > max_size:
            labels[labels == prop.label] = 0

    # Step 8: Fill holes
    # "After both thresholding and declumping"
    # Since there's no declumping (method=None), this just fills holes after thresholding
    for label_id in np.unique(labels):
        if label_id == 0:
            continue
        mask = labels == label_id
        filled = ndimage.binary_fill_holes(mask)
        labels[filled] = label_id

    # Step 9: Remove border objects
    labels = segmentation.clear_border(labels)

    # Step 10: Relabel consecutively
    labels = morphology.label(labels > 0)

    print(f"Detected {labels.max()} edge spots")

    return labels


# Test with synthetic data
if __name__ == "__main__":
    print("="*80)
    print("ROBUST BACKGROUND THRESHOLD ANALYSIS")
    print("="*80)

    # Create test image
    np.random.seed(42)
    test_image = np.random.rand(100, 100) * 0.1  # Background

    # Add some bright spots
    test_image[20:25, 20:25] = 0.8
    test_image[60:68, 60:68] = 0.6
    test_image[80:84, 30:34] = 0.9

    print("\nTest image statistics:")
    print(f"  Mean: {test_image.mean():.4f}")
    print(f"  Median: {np.median(test_image):.4f}")
    print(f"  Std: {test_image.std():.4f}")
    print(f"  Max: {test_image.max():.4f}")

    # Calculate threshold with Module 11 parameters
    print("\nCalculating threshold with Module 11 parameters:")
    print("  Lower outlier fraction: 0.3")
    print("  Upper outlier fraction: 0.1")
    print("  Averaging method: median")
    print("  Variance method: std")
    print("  Number of deviations: 6")
    print("  Correction factor: 2.0")
    print("  Lower bound: 0.0035")
    print("  Upper bound: 1.0")

    threshold = robust_background_threshold(
        test_image,
        lower_outlier_fraction=0.3,
        upper_outlier_fraction=0.1,
        averaging_method='median',
        variance_method='std',
        num_deviations=6,
        correction_factor=2.0
    )

    print(f"\nCalculated threshold: {threshold:.6f}")
    print(f"After bounds [0.0035, 1.0]: {np.clip(threshold, 0.0035, 1.0):.6f}")

    # Show what gets detected
    spots = identify_edge_spots_correct(
        test_image,
        min_diameter=5,
        max_diameter=80,
        threshold_smoothing_scale=1.3,
        threshold_correction_factor=2.0,
        lower_outlier_fraction=0.3,
        upper_outlier_fraction=0.1,
        num_deviations=6,
        lower_bound=0.0035,
        upper_bound=1.0
    )

    from skimage import measure as sk_measure
    print("\nSpot areas:")
    props = sk_measure.regionprops(spots)
    for i, prop in enumerate(props, 1):
        print(f"  Spot {i}: {prop.area} pixels, centroid at {prop.centroid}")

