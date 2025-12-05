# Module 11: IdentifyPrimaryObjects (edge_spots) - DEFINITIVE REFERENCE

## Overview

**Purpose:** Detect bright mitochondrial spots in the peripheral (non-perinuclear) region  
**Input:** `Subtract_perinucleus_Miro160mer` (MIRO160mer image with 15px perinuclear region masked out)  
**Output:** `edge_spots` (labeled image of detected spots)  
**Status:** ✅ ENABLED

---

## Complete Parameter List (Verified from HDF5)

### Basic Parameters
```
Input image: Subtract_perinucleus_Miro160mer
Output objects: edge_spots
Typical diameter (Min,Max): 5,80 pixels
Discard objects outside diameter range?: Yes
Discard border objects?: Yes
```

### Declumping Parameters
```
Method to distinguish clumped objects: None
Method to draw dividing lines: Intensity (not used since method=None)
Size of smoothing filter: 10 (not used since method=None)
Minimum distance between local maxima: 7.0 pixels (not used since method=None)
Speed up by using lower-resolution image?: Yes (not used since method=None)
Fill holes: After both thresholding and declumping
Auto-calculate smoothing filter?: Yes
Auto-calculate minimum distance?: Yes
Maximum objects: 500
Handling if excessive objects: Continue
```

**CRITICAL:** `Method to distinguish clumped objects: None` means:
- NO watershed segmentation
- NO declumping of touching objects
- Objects are identified purely by connected components
- All "declumping" parameters (smoothing filter, min distance) are IGNORED

### Advanced Thresholding Parameters
```
Use advanced settings?: Yes
Threshold setting version: 12
Threshold strategy: Global
Thresholding method: Robust Background
```

#### Robust Background Parameters
```
Threshold smoothing scale: 1.3
Threshold correction factor: 2
Lower threshold bound: 0.0035
Upper threshold bound: 1.0
Manual threshold: 0.05 (not used with Robust Background)
Select measurement to threshold with: None
```

#### Two-Class Thresholding
```
Two-class or three-class?: Two classes
Middle class assignment: Foreground (not applicable for two-class)
Log transform before thresholding?: No
```

#### Robust Background Outlier Parameters
```
Size of adaptive window: 50 (for adaptive methods - not used with Global)
Lower outlier fraction: 0.3
Upper outlier fraction: 0.1
Averaging method: Median
Variance method: Standard deviation
Number of deviations: 6
```

#### Fallback Method
```
Thresholding method: Minimum Cross-Entropy
```
(Used if Robust Background fails - rare)

---

## Robust Background Algorithm - EXACT IMPLEMENTATION

### Step-by-Step Process

**Step 1: Extract non-zero pixels**
```python
pixels = image.flatten()
pixels = pixels[pixels > 0]  # Ignore zero (masked) pixels
pixels_sorted = np.sort(pixels)
n_pixels = len(pixels_sorted)
```

**Step 2: Remove outliers**
```python
# Remove lowest 30% and highest 10% of pixel values
lower_idx = int(n_pixels * 0.3)  # Lower outlier fraction
upper_idx = int(n_pixels * 0.9)  # 1.0 - upper outlier fraction

robust_pixels = pixels_sorted[lower_idx:upper_idx]
```

**Why remove outliers?**
- Lower 30%: Removes background noise, masked regions
- Upper 10%: Removes very bright spots (the objects we want to detect)
- Remaining 60%: Represents the "robust background"

**Step 3: Calculate robust statistics**
```python
# Use MEDIAN (not mean) for robustness
robust_mean = np.median(robust_pixels)

# Use standard deviation of the robust pixels
robust_std = np.std(robust_pixels, ddof=1)  # Sample std
```

**Step 4: Calculate threshold**
```python
# Base threshold: mean + 6 standard deviations
threshold = robust_mean + (6 * robust_std)

# Apply correction factor
threshold *= 2.0

# Apply bounds
threshold = np.clip(threshold, 0.0035, 1.0)
```

**Why 6 standard deviations?**
- Very stringent (99.9999998% of normal distribution is below 6σ)
- Only detects VERY bright spots, far above background
- Then multiplied by 2 for even more stringency

**Step 5: Apply threshold**
```python
binary = image > threshold
```

**Step 6: Smooth the BINARY image**
```python
# CRITICAL: Smoothing is applied to the BINARY mask, not the grayscale image
binary_float = binary.astype(float)
binary_smoothed = filters.gaussian(binary_float, sigma=1.3)
binary = binary_smoothed > 0.5
```

**Why smooth the binary?**
- Threshold smoothing scale = 1.3
- Smooths jagged edges from thresholding
- Makes object boundaries smoother
- Different from pre-smoothing the image

---

## Complete Python Implementation

```python
import numpy as np
from skimage import filters, morphology, measure, segmentation
from scipy import ndimage

def robust_background_threshold(image, 
                                 lower_outlier_fraction=0.3,
                                 upper_outlier_fraction=0.1,
                                 averaging_method='median',
                                 variance_method='std',
                                 num_deviations=6,
                                 correction_factor=2.0,
                                 lower_bound=0.0035,
                                 upper_bound=1.0):
    """
    CellProfiler's Robust Background thresholding
    
    This is the EXACT algorithm used by CellProfiler.
    """
    # Extract non-zero pixels
    pixels = image.flatten()
    pixels = pixels[pixels > 0]
    
    if len(pixels) == 0:
        return 0.0
    
    # Sort pixels
    pixels_sorted = np.sort(pixels)
    n_pixels = len(pixels_sorted)
    
    # Calculate outlier indices
    lower_idx = int(n_pixels * lower_outlier_fraction)
    upper_idx = int(n_pixels * (1.0 - upper_outlier_fraction))
    
    # Get robust pixels (exclude outliers)
    robust_pixels = pixels_sorted[lower_idx:upper_idx]
    
    if len(robust_pixels) == 0:
        return lower_bound
    
    # Calculate robust mean (using median)
    if averaging_method.lower() == 'median':
        robust_mean = np.median(robust_pixels)
    else:
        robust_mean = np.mean(robust_pixels)
    
    # Calculate robust variance (using std)
    if variance_method.lower() == 'std':
        robust_std = np.std(robust_pixels, ddof=1)
    else:  # MAD
        deviations = np.abs(robust_pixels - np.median(robust_pixels))
        robust_std = np.median(deviations) * 1.4826
    
    # Calculate threshold
    threshold = robust_mean + (num_deviations * robust_std)
    
    # Apply correction factor
    threshold *= correction_factor
    
    # Apply bounds
    threshold = np.clip(threshold, lower_bound, upper_bound)
    
    return threshold


def identify_edge_spots(masked_image):
    """
    Module 11: Identify edge spots
    
    EXACT implementation matching CellProfiler parameters
    """
    img = masked_image.copy()
    
    # Parameters from Module 11
    MIN_DIAMETER = 5
    MAX_DIAMETER = 80
    THRESHOLD_SMOOTHING = 1.3
    CORRECTION_FACTOR = 2.0
    LOWER_OUTLIER = 0.3
    UPPER_OUTLIER = 0.1
    NUM_DEVIATIONS = 6
    LOWER_BOUND = 0.0035
    UPPER_BOUND = 1.0
    
    # Step 1: Calculate Robust Background threshold
    threshold = robust_background_threshold(
        img,
        lower_outlier_fraction=LOWER_OUTLIER,
        upper_outlier_fraction=UPPER_OUTLIER,
        averaging_method='median',
        variance_method='std',
        num_deviations=NUM_DEVIATIONS,
        correction_factor=CORRECTION_FACTOR,
        lower_bound=LOWER_BOUND,
        upper_bound=UPPER_BOUND
    )
    
    # Step 2: Apply threshold
    binary = img > threshold
    
    # Step 3: Smooth the BINARY image (threshold smoothing)
    binary_float = binary.astype(float)
    binary_smoothed = filters.gaussian(binary_float, sigma=THRESHOLD_SMOOTHING)
    binary = binary_smoothed > 0.5
    
    # Step 4: Size filtering (minimum)
    min_size = int(np.pi * (MIN_DIAMETER/2)**2)
    binary = morphology.remove_small_objects(binary, min_size=min_size)
    
    # Step 5: Label objects
    # CRITICAL: Method = None means NO watershed, just connected components
    labels = morphology.label(binary)
    
    # Step 6: Filter by maximum size
    max_size = int(np.pi * (MAX_DIAMETER/2)**2)
    props = measure.regionprops(labels)
    for prop in props:
        if prop.area > max_size:
            labels[labels == prop.label] = 0
    
    # Step 7: Fill holes
    # "After both thresholding and declumping"
    # Since there's no declumping (method=None), this is after thresholding only
    for label_id in np.unique(labels):
        if label_id == 0:
            continue
        mask = labels == label_id
        filled = ndimage.binary_fill_holes(mask)
        labels[filled] = label_id
    
    # Step 8: Remove border objects
    labels = segmentation.clear_border(labels)
    
    # Step 9: Relabel consecutively
    labels = morphology.label(labels > 0)
    
    return labels
```

---

## Key Insights

### 1. No Declumping
**Method to distinguish clumped objects: None**

This means:
- Touching spots are NOT separated
- Only simple connected component labeling
- Each connected region = one object
- Much simpler than Module 7 (nuclei) which uses watershed

### 2. Very Stringent Thresholding
The combination of:
- 6 standard deviations (extreme outlier)
- Correction factor of 2.0 (doubles threshold)
- Removing lower 30% and upper 10% of pixels

Results in: **ONLY the brightest spots are detected**

### 3. Robust to Background Variation
By using:
- Median instead of mean (robust to outliers)
- Outlier fraction removal (ignores extremes)
- Adaptive to image content (not fixed threshold)

### 4. Binary Smoothing (Not Image Smoothing)
**Threshold smoothing scale: 1.3**

- Applied AFTER thresholding
- Smooths the binary mask, not the image
- Creates smoother object boundaries
- Different from Gaussian pre-filtering

---

## Comparison: Module 7 (Nuclei) vs Module 11 (edge_spots)

| Feature | Module 7 (Nuclei) | Module 11 (edge_spots) |
|---------|------------------|------------------------|
| **Input** | Hoechst | Subtract_perinucleus_Miro160mer |
| **Diameter** | 30-100 px | 5-80 px |
| **Thresholding** | 3-class Otsu | Robust Background |
| **Correction** | 0.45 | 2.0 |
| **Declumping** | Intensity (watershed) | None (connected components) |
| **Smoothing** | Pre-threshold (σ=2) | Post-threshold (σ=1.3) |
| **Stringency** | Moderate (captures most nuclei) | Very high (only brightest spots) |
| **Purpose** | Segment all nuclei | Detect bright peripheral mitochondria |

---

## Troubleshooting

### Too many spots detected
**Adjust:**
1. Increase `Threshold correction factor` (currently 2.0) → 2.5 or 3.0
2. Increase `Number of deviations` (currently 6) → 7 or 8
3. Increase `Lower bound` (currently 0.0035) → 0.005 or 0.01
4. Increase `Lower outlier fraction` (currently 0.3) → 0.4

### Too few spots detected
**Adjust:**
1. Decrease `Threshold correction factor` (currently 2.0) → 1.5 or 1.0
2. Decrease `Number of deviations` (currently 6) → 4 or 5
3. Decrease `Lower bound` (currently 0.0035) → 0.002 or 0.001
4. Decrease `Lower outlier fraction` (currently 0.3) → 0.2

### Spots have jagged edges
**Adjust:**
1. Increase `Threshold smoothing scale` (currently 1.3) → 2.0 or 2.5

### Need to separate touching spots
**Change:**
- `Method to distinguish clumped objects` from "None" to "Intensity"
- This will enable watershed declumping (not currently used)

---

## Mathematical Formula Summary

```
STEP 1: Sort non-zero pixels
pixels_sorted = sort(pixels[pixels > 0])

STEP 2: Extract robust pixels
n = len(pixels_sorted)
robust_pixels = pixels_sorted[int(n*0.3) : int(n*0.9)]

STEP 3: Calculate robust statistics
robust_mean = median(robust_pixels)
robust_std = std(robust_pixels)

STEP 4: Calculate threshold
threshold = robust_mean + (6 * robust_std)
threshold = threshold * 2.0
threshold = clip(threshold, 0.0035, 1.0)

STEP 5: Apply threshold and smooth
binary = image > threshold
binary = gaussian(binary, σ=1.3) > 0.5

STEP 6: Size filter and clean up
binary = remove_small_objects(binary, min_size)
labels = connected_components(binary)
labels = remove_large_objects(labels, max_size)
labels = fill_holes(labels)
labels = clear_border(labels)
```

---

## Test Case

```python
# Example calculation for typical image
# Assuming background pixels ~0.05, bright spots ~0.8

# Pixel distribution after masking:
# Background: 0.03 - 0.08 (90% of pixels)
# Bright spots: 0.5 - 0.9 (5% of pixels)
# Masked (zero): 0.0 (5% of pixels)

# Step 1: Sort and remove outliers
# After removing lower 30% and upper 10%:
# Robust pixels: 0.04 - 0.07 (60% of original)

# Step 2: Calculate statistics
robust_mean = median([0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07])
            ≈ 0.055

robust_std = std([0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07])
           ≈ 0.011

# Step 3: Calculate threshold
threshold = 0.055 + (6 * 0.011)
          = 0.055 + 0.066
          = 0.121

threshold *= 2.0
          = 0.242

# Step 4: Apply bounds
threshold = clip(0.242, 0.0035, 1.0)
          = 0.242

# Result: Only pixels > 0.242 are detected
# This captures the bright spots (0.5-0.9) but not background (0.03-0.08)
```

---

## References

- CellProfiler Manual: https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.4/modules/imageprocessing.html#identifyprimaryobjects
- Source Code: https://github.com/CellProfiler/CellProfiler/blob/master/src/subpackages/core/cellprofiler_core/threshold/_robust_background_threshold.py
- Paper: Jones et al., "CellProfiler Analyst: data exploration and analysis software for complex image-based screens", BMC Bioinformatics 2008

---

**Last verified:** Against HDF5 and JSON exports from `231115_combined_pipeline_new_nomito_fixed_moments.cpproj`  
**CellProfiler version:** 4.2.4  
**Module revision:** 15
