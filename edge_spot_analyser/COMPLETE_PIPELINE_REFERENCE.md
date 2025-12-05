# CellProfiler Pipeline - Complete Technical Reference

**Pipeline File:** `231115_combined_pipeline_new_nomito_fixed_moments.cpproj`  
**CellProfiler Version:** 4.2.4  
**Pipeline Version:** v6  
**Date Revision:** 424  
**Total Modules:** 22 (17 enabled)

---

## Pipeline Purpose

Analyze mitochondrial (MIRO160mer) distribution around nuclei, with specific focus on:
1. Nuclear segmentation from Hoechst staining
2. Perinuclear region definition
3. Peripheral mitochondrial spot detection
4. Statistical characterization of mitochondrial distribution

---

## Input Configuration

### Module 1: Images
**Status:** ✅ Enabled  
**Purpose:** File import configuration

```
Filter images?: Images only
Rule criteria: and (extension does isimage) (directory doesnot startwith ".")
```

**Python equivalent:**
```python
from pathlib import Path
import imageio

def load_image_files(directory):
    valid_extensions = {'.tif', '.tiff', '.png', '.jpg', '.jpeg'}
    files = [f for f in Path(directory).rglob('*') 
             if f.suffix.lower() in valid_extensions 
             and not f.parent.name.startswith('.')]
    return files
```

---

### Module 2: Metadata
**Status:** ✅ Enabled (but extraction disabled)  
**Purpose:** Extract metadata from filenames

```
Extract metadata?: No
Metadata extraction method: Extract from file/folder names
Regular expression (filename): ^(?P<Type>.+)_(?P<Construct>\w+).*$
Regular expression (folder): ^(?P<Type>.*)_(?P<Construct>\w+)
Filter criteria: and (file does contain "TRAK2700")
```

**Note:** Metadata extraction is disabled, so this step doesn't affect the pipeline.

---

### Module 3: NamesAndTypes
**Status:** ✅ Enabled  
**Purpose:** Assign channel names based on filename patterns

```
Assign a name to: Images matching rules
Image set matching method: Order
Maximum intensity: 255.0
Process as 3D?: No
Pixel spacing (X,Y,Z): 1.0, 1.0, 1.0

Channel 1 - Hoechst (Nuclear):
  File pattern: file ends with "405.tif"
  Name: Hoechst
  Type: Grayscale image
  Intensity range: From image metadata

Channel 2 - MIRO160mer (Mitochondrial):
  File pattern: file ends with "488.tif"
  Name: MIRO160mer
  Type: Grayscale image
  Intensity range: From image metadata
```

**Python equivalent:**
```python
def assign_channels(file_list):
    channels = {'Hoechst': [], 'MIRO160mer': []}
    for file in file_list:
        if file.name.endswith('405.tif'):
            channels['Hoechst'].append(file)
        elif file.name.endswith('488.tif'):
            channels['MIRO160mer'].append(file)
    return channels
```

---

### Module 4: Groups
**Status:** ✅ Enabled (but grouping disabled)  
**Purpose:** Batch processing configuration

```
Do you want to group your images?: No
```

**Note:** No grouping is used, so images are processed independently.

---

## Preprocessing (DISABLED)

### Module 5: EnhanceOrSuppressFeatures
**Status:** ❌ DISABLED  
Would enhance speckle features in Hoechst image if enabled.

### Module 6: RescaleIntensity
**Status:** ❌ DISABLED  
Would rescale MIRO160mer intensity if enabled.

---

## SEGMENTATION PIPELINE

### Module 7: IdentifyPrimaryObjects → **Nuclei**
**Status:** ✅ Enabled  
**Purpose:** Segment nuclei from Hoechst channel

#### Basic Parameters
```
Input image: Hoechst
Output objects: Nuclei
Typical diameter: 30-100 pixels
Discard objects outside diameter range?: Yes
Discard border objects?: Yes
```

#### Declumping Parameters
```
Method to distinguish clumped objects: Intensity
Method to draw dividing lines: Intensity
Size of smoothing filter: 10
Minimum distance between local maxima: 30 pixels
Speed up by using lower-resolution image?: Yes
Fill holes: After declumping only
Auto-calculate smoothing filter?: Yes
Auto-calculate minimum distance?: No
Maximum number of objects: 500
Handling if excessive objects: Continue
```

#### Thresholding Parameters (CRITICAL)
```
Use advanced settings?: Yes
Threshold setting version: 12
Threshold strategy: Global
Thresholding method: Otsu

IMPORTANT - Three-class thresholding:
  Two-class or three-class?: Three classes
  Middle class assignment: Background
  
Threshold smoothing scale: 2
Threshold correction factor: 0.45
Lower threshold bound: 0.00195
Upper threshold bound: 0.01
Manual threshold: 0.0
Log transform before thresholding?: No

Outlier handling:
  Lower outlier fraction: 0.05
  Upper outlier fraction: 0.05
  Averaging method: Mean
  Variance method: Standard deviation
  Number of deviations: 2
```

#### Python Implementation
```python
from skimage import filters, morphology, measure, segmentation
from scipy import ndimage
import numpy as np

def identify_nuclei(hoechst_img):
    """
    Segment nuclei using three-class Otsu thresholding
    """
    # 1. Smooth image (threshold smoothing scale = 2)
    smoothed = filters.gaussian(hoechst_img, sigma=2)
    
    # 2. Three-class Otsu thresholding
    from skimage.filters import threshold_multiotsu
    thresholds = threshold_multiotsu(smoothed, classes=3)
    
    # Use higher threshold and apply correction factor
    thresh_value = thresholds[1] * 0.45
    
    # 3. Apply threshold bounds
    thresh_value = np.clip(thresh_value, 0.00195, 0.01)
    
    # 4. Create binary mask
    binary = smoothed > thresh_value
    
    # 5. Remove small/large objects
    min_size = int(np.pi * (30/2)**2)  # ~706 pixels
    max_size = int(np.pi * (100/2)**2)  # ~7854 pixels
    binary = morphology.remove_small_objects(binary, min_size=min_size)
    
    # 6. Watershed segmentation for declumping
    # Distance transform
    distance = ndimage.distance_transform_edt(binary)
    
    # Find local maxima (min distance = 30 pixels)
    # Smooth distance with sigma=10 for smoothing filter
    distance_smooth = filters.gaussian(distance, sigma=10)
    local_max = morphology.local_maxima(distance_smooth)
    
    # Create markers from local maxima, enforcing min distance
    coords = np.argwhere(local_max)
    markers = np.zeros_like(binary, dtype=int)
    label_id = 1
    for i, coord in enumerate(coords):
        # Check minimum distance to other markers
        if markers[tuple(coord)] == 0:
            y, x = coord
            # Mark region
            markers[y, x] = label_id
            # Suppress nearby maxima (within 30 pixels)
            yy, xx = np.ogrid[-30:31, -30:31]
            circle = (yy**2 + xx**2 <= 30**2)
            y_min = max(0, y-30)
            y_max = min(markers.shape[0], y+31)
            x_min = max(0, x-30)
            x_max = min(markers.shape[1], x+31)
            circle_crop = circle[
                (30-(y-y_min)):(30+(y_max-y-1)),
                (30-(x-x_min)):(30+(x_max-x-1))
            ]
            markers[y_min:y_max, x_min:x_max][circle_crop] = label_id
            label_id += 1
    
    markers = morphology.label(markers > 0)
    
    # 7. Watershed
    labels = segmentation.watershed(-distance, markers, mask=binary)
    
    # 8. Filter by size
    for region in measure.regionprops(labels):
        if region.area < min_size or region.area > max_size:
            labels[labels == region.label] = 0
    
    # 9. Fill holes (after declumping only)
    for label_id in np.unique(labels):
        if label_id == 0:
            continue
        mask = labels == label_id
        filled = ndimage.binary_fill_holes(mask)
        labels[filled] = label_id
    
    # 10. Remove border objects
    labels = segmentation.clear_border(labels)
    
    # 11. Relabel consecutively
    labels = morphology.label(labels > 0)
    
    return labels
```

---

### Module 8: ExpandOrShrinkObjects → **Expand_Nuclei**
**Status:** ✅ Enabled  
**Purpose:** Create perinuclear region for measurements

```
Input objects: Nuclei
Output objects: Expand_Nuclei
Operation: Expand objects by a specified number of pixels
Number of pixels: 10
Fill holes to single point?: No
```

#### Python Implementation
```python
def expand_nuclei(nuclei_labels, expansion_pixels=10):
    """
    Expand each nucleus by specified number of pixels
    Preserves individual labels
    """
    from skimage.segmentation import expand_labels
    expanded = expand_labels(nuclei_labels, distance=expansion_pixels)
    return expanded
```

---

### Module 9: ExpandOrShrinkObjects → **Expand_Nuclei_for_mask**
**Status:** ✅ Enabled  
**Purpose:** Create larger region for masking perinuclear area

```
Input objects: Nuclei
Output objects: Expand_Nuclei_for_mask
Operation: Expand objects by a specified number of pixels
Number of pixels: 15
Fill holes to single point?: No
```

#### Python Implementation
```python
def expand_nuclei_for_mask(nuclei_labels, expansion_pixels=15):
    """
    Expand each nucleus by 15 pixels for masking
    """
    from skimage.segmentation import expand_labels
    expanded = expand_labels(nuclei_labels, distance=expansion_pixels)
    return expanded
```

---

### Module 10: MaskImage → **Subtract_perinucleus_Miro160mer**
**Status:** ✅ Enabled  
**Purpose:** Remove perinuclear region to isolate peripheral mitochondria

```
Input image: MIRO160mer
Output image: Subtract_perinucleus_Miro160mer
Use objects or image as mask?: Objects
Mask object: Expand_Nuclei_for_mask
Invert the mask?: Yes
```

**Effect:** Sets pixels within 15px of nuclei to 0, keeping only peripheral regions

#### Python Implementation
```python
def mask_perinuclear_region(miro_image, expand_nuclei_for_mask):
    """
    Mask out the perinuclear region (inverted mask)
    Keep only pixels OUTSIDE the expanded nuclear region
    """
    mask = expand_nuclei_for_mask > 0
    masked = miro_image.copy()
    masked[mask] = 0  # Inverted: remove perinuclear, keep periphery
    return masked
```

---

### Module 11: IdentifyPrimaryObjects → **edge_spots**
**Status:** ✅ Enabled  
**Purpose:** Detect bright mitochondrial spots in peripheral regions

#### Basic Parameters
```
Input image: Subtract_perinucleus_Miro160mer
Output objects: edge_spots
Typical diameter: 5-80 pixels
Discard objects outside diameter range?: Yes
Discard border objects?: Yes
```

#### Declumping Parameters
```
Method to distinguish clumped objects: None
Method to draw dividing lines: Intensity
Size of smoothing filter: 10
Minimum distance between local maxima: 7.0 pixels
Speed up by using lower-resolution image?: Yes
Fill holes: After both thresholding and declumping
Auto-calculate smoothing filter?: Yes
Auto-calculate minimum distance?: Yes
Maximum number of objects: 500
Handling if excessive objects: Continue
```

#### Thresholding Parameters (CRITICAL)
```
Use advanced settings?: Yes
Threshold setting version: 12
Threshold strategy: Global
Thresholding method: Robust Background

ROBUST BACKGROUND SPECIFIC PARAMETERS:
  Threshold smoothing scale: 1.3
  Threshold correction factor: 2
  Lower threshold bound: 0.0035
  Upper threshold bound: 1.0
  Manual threshold: 0.05
  
  Two-class or three-class?: Two classes
  Middle class assignment: Foreground
  Log transform before thresholding?: No
  
  Adaptive window parameters:
    Size of adaptive window: 50
    Lower outlier fraction: 0.3
    Upper outlier fraction: 0.1
    Averaging method: Median
    Variance method: Standard deviation
    Number of deviations: 6
    
  Fallback thresholding method: Minimum Cross-Entropy
```

#### Robust Background Algorithm Explained

The Robust Background thresholding is a CellProfiler-specific method:

1. **Calculate robust statistics using outlier fractions:**
   - Remove lower 30% and upper 10% of non-zero pixels
   - Use median of remaining pixels as robust mean
   - Calculate robust standard deviation

2. **Compute threshold:**
   ```
   threshold = robust_mean + (num_deviations * robust_std)
   threshold *= correction_factor
   threshold = max(threshold, lower_bound)
   threshold = min(threshold, upper_bound)
   ```

3. **Apply smoothing to binary result** with sigma=1.3

#### Python Implementation
```python
def identify_edge_spots(masked_miro_image):
    """
    Detect bright spots using Robust Background thresholding
    """
    img = masked_miro_image
    
    # Get non-zero pixels only
    img_nonzero = img[img > 0]
    
    if len(img_nonzero) == 0:
        return np.zeros_like(img, dtype=int)
    
    # Robust Background thresholding
    # 1. Remove outliers (lower 30%, upper 10%)
    lower_percentile = np.percentile(img_nonzero, 30)
    upper_percentile = np.percentile(img_nonzero, 90)
    
    # Keep only middle 60%
    robust_pixels = img_nonzero[
        (img_nonzero >= lower_percentile) & 
        (img_nonzero <= upper_percentile)
    ]
    
    # 2. Calculate robust statistics
    robust_mean = np.median(robust_pixels)  # Median as robust mean
    robust_std = np.std(robust_pixels)
    
    # 3. Calculate threshold with 6 deviations
    threshold = robust_mean + (6 * robust_std)
    
    # 4. Apply correction factor
    threshold *= 2.0
    
    # 5. Apply bounds
    threshold = max(threshold, 0.0035)
    threshold = min(threshold, 1.0)
    
    # 6. Create binary mask
    binary = img > threshold
    
    # 7. Smooth binary result (threshold smoothing scale = 1.3)
    binary_smooth = filters.gaussian(binary.astype(float), sigma=1.3)
    binary = binary_smooth > 0.5
    
    # 8. Size filtering
    min_size = int(np.pi * (5/2)**2)   # ~20 pixels
    max_size = int(np.pi * (80/2)**2)  # ~5027 pixels
    binary = morphology.remove_small_objects(binary, min_size=min_size)
    
    # 9. Label objects (no declumping - method is "None")
    labels = morphology.label(binary)
    
    # 10. Filter by size
    for region in measure.regionprops(labels):
        if region.area > max_size:
            labels[labels == region.label] = 0
    
    # 11. Fill holes
    for label_id in np.unique(labels):
        if label_id == 0:
            continue
        mask = labels == label_id
        filled = ndimage.binary_fill_holes(mask)
        labels[filled] = label_id
    
    # 12. Remove border objects
    labels = segmentation.clear_border(labels)
    
    # 13. Relabel
    labels = morphology.label(labels > 0)
    
    return labels
```

**Key Differences from Simple Thresholding:**
- Robust to outliers (ignores extreme values)
- Adaptive to local background
- Uses median instead of mean for robustness
- Large number of standard deviations (6) means only very bright spots
- Correction factor of 2 further increases stringency

---

## MEASUREMENT MODULES

### Module 12: MeasureObjectIntensity
**Status:** ✅ Enabled  
**Purpose:** Measure intensity of edge spots

```
Images to measure: Subtract_perinucleus_Miro160mer
Objects to measure: edge_spots
```

**Measurements computed:**
- Mean intensity
- Median intensity
- Std intensity
- Min/Max intensity
- Integrated intensity (sum)
- Lower/Upper quartile
- MAD (Median Absolute Deviation)
- Mean/Std intensity at edge
- Mass displacement

---

### Module 13: FilterObjects → **FilterObjects**
**Status:** ✅ Enabled  
**Purpose:** Filter spots by edge intensity quality

```
Input objects: edge_spots
Output objects: FilterObjects
Filtering mode: Measurements
Filtering method: Limits
Measurement: Intensity_MeanIntensityEdge_Subtract_perinucleus_Miro160mer
Minimum value: 0.0
Maximum value: 1.0
Keep removed objects?: No
```

**Effect:** Removes spots with edge intensity outside [0.0, 1.0] range (essentially keeps all spots with valid intensities)

#### Python Implementation
```python
def filter_edge_spots(edge_spots_labels, edge_intensities_df):
    """
    Filter spots by mean edge intensity
    """
    filtered = np.zeros_like(edge_spots_labels)
    
    for _, row in edge_intensities_df.iterrows():
        label_id = int(row['label'])
        mean_edge = row['mean_intensity_edge']
        
        # Keep if within bounds
        if 0.0 <= mean_edge <= 1.0:
            filtered[edge_spots_labels == label_id] = label_id
    
    # Relabel consecutively
    filtered = morphology.label(filtered > 0)
    return filtered
```

---

### Module 14: IdentifyTertiaryObjects → **Perinuclear_region**
**Status:** ✅ Enabled  
**Purpose:** Create annular perinuclear region

```
Larger objects: Expand_Nuclei
Smaller objects: Nuclei
Output objects: Perinuclear_region
Shrink smaller object prior to subtraction?: Yes
```

**Effect:** Creates a ring-shaped region between the nucleus and the expanded nucleus (10-pixel annulus around each nucleus)

#### Python Implementation
```python
def create_perinuclear_region(nuclei_labels, expand_nuclei_labels):
    """
    Create perinuclear ring by subtracting nucleus from expanded nucleus
    Shrinks nucleus by 1 pixel before subtraction
    """
    # Shrink nuclei slightly
    shrunk_nuclei = morphology.erosion(nuclei_labels > 0, morphology.disk(1))
    
    # Create ring
    perinuclear = expand_nuclei_labels.copy()
    for label_id in np.unique(expand_nuclei_labels):
        if label_id == 0:
            continue
        # Remove shrunk nucleus from this expanded region
        mask = expand_nuclei_labels == label_id
        mask[shrunk_nuclei & (nuclei_labels == label_id)] = False
        perinuclear[~mask & (expand_nuclei_labels == label_id)] = 0
    
    return perinuclear
```

---

### Module 15: MeasureObjectIntensity
**Status:** ✅ Enabled  
**Purpose:** Measure nuclear intensity

```
Images to measure: Hoechst
Objects to measure: Nuclei
```

**Measurements:** All standard intensity metrics for Hoechst channel in nuclear regions

---

### Module 16: MeasureObjectIntensity
**Status:** ✅ Enabled  
**Purpose:** Measure mitochondrial intensity in nuclear regions

```
Images to measure: MIRO160mer
Objects to measure: Expand_Nuclei, Perinuclear_region
```

**Measurements:** All standard intensity metrics for MIRO160mer channel in:
- Expanded nuclear regions
- Perinuclear rings

---

### Module 17: CalculateMoments
**Status:** ✅ Enabled  
**Purpose:** Calculate distribution statistics

```
Image to measure: MIRO160mer
Objects to measure: 
  - Perinuclear_region
  - Expand_Nuclei
Moments to compute: Mean, Standard Deviation, Skewness, Kurtosis
```

#### Moment Definitions

**Mean (μ):**
```
μ = Σ(xi) / n
```
Average intensity value

**Standard Deviation (σ):**
```
σ = sqrt(Σ(xi - μ)² / n)
```
Spread of intensity distribution

**Skewness:**
```
skewness = E[(X - μ)³] / σ³
```
- Positive: tail extends to higher intensities (bright spots)
- Negative: tail extends to lower intensities
- ~0: symmetric distribution

**Kurtosis:**
```
kurtosis = E[(X - μ)⁴] / σ⁴ - 3
```
- Positive: heavy tails, sharp peak (concentrated intensity)
- Negative: light tails, flat peak (uniform intensity)
- ~0: normal distribution

#### Python Implementation
```python
from scipy import stats

def calculate_moments(image, labels):
    """
    Calculate statistical moments for each labeled region
    """
    results = []
    
    for region in measure.regionprops(labels, intensity_image=image):
        intensities = image[labels == region.label]
        
        results.append({
            'label': region.label,
            'mean': np.mean(intensities),
            'std': np.std(intensities, ddof=1),  # Sample std
            'skewness': stats.skew(intensities),
            'kurtosis': stats.kurtosis(intensities, fisher=True)  # Excess kurtosis
        })
    
    return pd.DataFrame(results)
```

---

### Module 18: CalculateHistogram
**Status:** ❌ DISABLED  
Would calculate 50-bin histograms if enabled

---

### Module 19: CalculateGini
**Status:** ✅ Enabled  
**Purpose:** Calculate Gini coefficient (inequality measure)

```
Image to measure: MIRO160mer
Objects to measure:
  - Expand_Nuclei
  - Perinuclear_region
```

#### Gini Coefficient Explained

The Gini coefficient measures inequality in intensity distribution:
- **0.0:** Perfectly uniform (all pixels same intensity)
- **1.0:** Perfectly unequal (all intensity in one pixel)

**Biological interpretation:**
- **Low Gini (0.0-0.3):** Evenly distributed mitochondria
- **High Gini (0.6-1.0):** Concentrated/clustered mitochondria

#### Formula
```
G = (2 * Σ(i * x_i)) / (n * Σ(x_i)) - (n + 1) / n

where:
  x_i = sorted intensity values
  i = rank (1 to n)
  n = number of pixels
```

#### Python Implementation
```python
def calculate_gini(image, labels):
    """
    Calculate Gini coefficient for intensity distribution
    """
    results = []
    
    for region in measure.regionprops(labels, intensity_image=image):
        intensities = image[labels == region.label]
        
        # Sort intensities
        sorted_intensities = np.sort(intensities)
        n = len(sorted_intensities)
        
        # Calculate cumulative sum
        cumsum = np.cumsum(sorted_intensities)
        
        # Gini coefficient
        if cumsum[-1] == 0:
            gini = 0.0
        else:
            gini = (2 * np.sum((np.arange(1, n+1)) * sorted_intensities)) / (n * cumsum[-1]) - (n + 1) / n
        
        results.append({
            'label': region.label,
            'gini': gini
        })
    
    return pd.DataFrame(results)
```

---

### Module 20: MeasureImageIntensity
**Status:** ✅ Enabled  
**Purpose:** Global image statistics

```
Images to measure: 
  - MIRO160mer
  - Subtract_perinucleus_Miro160mer
Measure only from enclosed objects?: No
Calculate custom percentiles?: No
```

**Measurements:**
- Mean, median, std, min, max
- Total intensity
- Total area
- Lower/Upper quartile
- MAD
- Percent maximal (% pixels at max intensity)

---

### Module 21: SaveImages
**Status:** ❌ DISABLED  
Would save MIRO160mer images if enabled

---

### Module 22: ExportToSpreadsheet
**Status:** ✅ Enabled  
**Purpose:** Export all measurements to CSV files

```
Column delimiter: Comma (",")
Add image metadata?: No
Add file and folder names?: Yes
Select measurements to export: Yes (specific measurements selected)
Output location: Z:\Christina\export_issue_test\moments_test\post_revision\241223_siRNA_motor_JNK_wt_DRH
Create GenePattern GCT file?: No
Representation of Nan/Inf: NaN
Add prefix to filenames?: No
Overwrite without warning?: No
```

#### Export Structure

**File 1: All_measurements.csv**
- Combines: Nuclei + edge_spots

**Files 2-6: Object-specific CSVs** (named after object type)
- edge_spots.csv
- Nuclei.csv
- Expand_Nuclei.csv
- Perinuclear_region.csv
- Image.csv

#### Measurements Exported

**For Nuclei:**
- Intensity_MedianIntensity_Hoechst
- Intensity_MeanIntensity_Hoechst
- Intensity_StdIntensity_Hoechst
- Intensity_IntegratedIntensity_Hoechst
- Number_Object_Number

**For Expand_Nuclei:**
- All intensity metrics (mean, median, std, min, max, quartiles, MAD, integrated)
- Edge intensity metrics
- Mass displacement
- All 50 histogram bins (even though module 18 is disabled)
- Moments (mean, std, skewness, kurtosis)
- Gini coefficient

**For Perinuclear_region:**
- All intensity metrics
- Edge intensity metrics  
- Mass displacement
- All 50 histogram bins
- Moments (mean, std, skewness, kurtosis)
- Gini coefficient

**For edge_spots:**
- All intensity metrics
- Edge intensity metrics
- Mass displacement
- Number_Object_Number

**For Image:**
- All global intensity statistics for both MIRO160mer images

---

## Complete Python Pipeline

```python
import numpy as np
from skimage import io, filters, measure, morphology, segmentation
from scipy import ndimage, stats
import pandas as pd

class CellProfilerPipelineComplete:
    
    def run_complete_pipeline(self, hoechst_path, miro_path):
        """
        Execute complete pipeline with exact CellProfiler parameters
        """
        # Load images
        hoechst = io.imread(hoechst_path, as_gray=True)
        miro = io.imread(miro_path, as_gray=True)
        
        # Normalize
        hoechst = hoechst.astype(float) / hoechst.max()
        miro = miro.astype(float) / miro.max()
        
        # Step 7: Segment nuclei
        nuclei = self.identify_nuclei(hoechst)
        
        # Steps 8-9: Expand nuclei
        expand_nuclei = self.expand_labels(nuclei, 10)
        expand_for_mask = self.expand_labels(nuclei, 15)
        
        # Step 10: Mask perinuclear region
        miro_masked = miro.copy()
        miro_masked[expand_for_mask > 0] = 0
        
        # Step 11: Identify edge spots
        edge_spots = self.identify_edge_spots(miro_masked)
        
        # Step 12: Measure edge spots
        edge_measurements = self.measure_intensity(edge_spots, miro_masked)
        
        # Step 13: Filter edge spots
        edge_spots_filtered = self.filter_spots(edge_spots, edge_measurements)
        
        # Step 14: Create perinuclear region
        perinuclear = self.create_perinuclear_region(nuclei, expand_nuclei)
        
        # Steps 15-16: Measure intensities
        nuclei_hoechst = self.measure_intensity(nuclei, hoechst)
        expand_miro = self.measure_intensity(expand_nuclei, miro)
        perinuclear_miro = self.measure_intensity(perinuclear, miro)
        
        # Step 17: Calculate moments
        expand_moments = self.calculate_moments(expand_nuclei, miro)
        perinuclear_moments = self.calculate_moments(perinuclear, miro)
        
        # Step 19: Calculate Gini
        expand_gini = self.calculate_gini(expand_nuclei, miro)
        perinuclear_gini = self.calculate_gini(perinuclear, miro)
        
        # Step 20: Image measurements
        image_stats_miro = self.measure_image(miro)
        image_stats_masked = self.measure_image(miro_masked)
        
        return {
            'objects': {
                'nuclei': nuclei,
                'expand_nuclei': expand_nuclei,
                'expand_for_mask': expand_for_mask,
                'edge_spots': edge_spots_filtered,
                'perinuclear': perinuclear
            },
            'measurements': {
                'nuclei_hoechst': nuclei_hoechst,
                'expand_miro': expand_miro,
                'perinuclear_miro': perinuclear_miro,
                'edge_spots': edge_measurements,
                'expand_moments': expand_moments,
                'perinuclear_moments': perinuclear_moments,
                'expand_gini': expand_gini,
                'perinuclear_gini': perinuclear_gini,
                'image_miro': image_stats_miro,
                'image_masked': image_stats_masked
            }
        }
    
    # [Include all the individual methods defined above]
```

---

## Critical Parameters Summary

| Module | Parameter | Value | Impact |
|--------|-----------|-------|--------|
| Nuclei segmentation | Diameter range | 30-100 px | Nucleus size filter |
| Nuclei segmentation | Thresholding | 3-class Otsu | Better separation |
| Nuclei segmentation | Correction factor | 0.45 | Lower = more nuclei |
| Nuclei segmentation | Min distance | 30 px | Prevents over-splitting |
| Nuclei expansion | Expansion | 10 px | Perinuclear ring width |
| Mask expansion | Expansion | 15 px | Excludes close mitochondria |
| Edge spots | Diameter range | 5-80 px | Spot size filter |
| Edge spots | Threshold method | Robust Background | Outlier-resistant |
| Edge spots | Correction factor | 2.0 | Higher = fewer spots |
| Edge spots | Lower bound | 0.0035 | Minimum intensity |
| Edge spots | Outlier fractions | 30%, 10% | Background estimation |
| Edge spots | Num deviations | 6 | Very stringent |

---

## Troubleshooting Guide

### Problem: Too many/few nuclei detected
**Adjust:**
- Threshold correction factor (0.45): Lower = more nuclei
- Diameter range (30-100): Adjust to your cell sizes
- Threshold bounds (0.00195-0.01): Adjust to your intensities

### Problem: No edge spots detected
**Adjust:**
- Threshold correction factor (2.0): Lower = more spots
- Lower bound (0.0035): Lower = dimmer spots detected
- Num deviations (6): Lower = more lenient threshold

### Problem: Nuclei are over-segmented
**Adjust:**
- Minimum distance between maxima (30 px): Increase
- Smoothing filter size (10): Increase

### Problem: Perinuclear ring too narrow/wide
**Adjust:**
- Expansion pixels (10): Change ring width
- Consider if 15px mask is appropriate for your data

---

## Output Files

When run in CellProfiler, this pipeline generates:

1. **All_measurements.csv** - Combined nuclei and edge_spots data
2. **edge_spots.csv** - Filtered spot measurements
3. **Nuclei.csv** - Nuclear measurements
4. **Expand_Nuclei.csv** - Expanded nuclear region measurements
5. **Perinuclear_region.csv** - Perinuclear ring measurements
6. **Image.csv** - Global image statistics

Each CSV contains measurements for all objects detected in all images processed.

---

## Version Information

```
CellProfiler: 4.2.4
Pipeline version: v6
Date revision: 424
Python libraries used:
  - scikit-image >= 0.19
  - scipy >= 1.7
  - numpy >= 1.21
  - pandas >= 1.3
```

---

## Notes

1. **Histogram bins**: Module 18 is disabled but histogram measurements are still exported - these will be empty or invalid
2. **Robust Background**: This is a CellProfiler-specific algorithm - the Python approximation may give slightly different results
3. **Edge intensity**: Calculated by dilating object masks and measuring intensity in the dilation-minus-object region
4. **Mass displacement**: Measures how far the intensity centroid is from the geometric centroid
5. **Three-class Otsu**: Assigns pixels to background, intermediate, or foreground classes - pipeline uses only foreground class

---

*This reference provides exact parameter values from the JSON export. For implementation questions, refer to the CellProfiler documentation at https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.4/index.html*
