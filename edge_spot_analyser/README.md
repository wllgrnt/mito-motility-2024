# Edge Spot Analyser

Python port of the CellProfiler pipeline for analyzing mitochondrial distribution around nuclei.

## Overview

This pipeline analyzes the distribution of mitochondria (MIRO160mer) around nuclei in microscopy images. It identifies bright peripheral mitochondrial spots and quantifies their spatial distribution using metrics like Gini coefficient, mass displacement, and coefficient of variation.

### What it does:

1. **Segments nuclei** from Hoechst (405nm) channel using 3-class Otsu thresholding
2. **Defines perinuclear regions** (10-15 pixel rings around nuclei)
3. **Detects bright peripheral mitochondrial spots** in MIRO160mer (488nm) channel
4. **Quantifies distribution metrics**: Gini coefficient, coefficient of variation (CoV), mass displacement, moments (skewness, kurtosis)
5. **Exports CellProfiler-compatible CSV files** for downstream analysis

## Installation

### Using uv (recommended)

```bash
# Sync dependencies and install the package
uv sync

# Or sync with optional dependencies
uv sync --extra yaml
uv sync --extra dev
uv sync --all-extras  # Install all optional dependencies
```

### Using pip (alternative)

```bash
pip install -e .
```

### Requirements

The package requires Python >=3.9 and will automatically install:
- numpy >= 1.20.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- scikit-image >= 0.18.0

## Input Data Format

The pipeline expects TIF files organized in date-based subdirectories:

```
inputs/
  ├── 241217/
  │   ├── WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_405.tif
  │   ├── WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_488.tif
  │   ├── WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY2_405.tif
  │   ├── WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY2_488.tif
  │   └── ...
  ├── 231120/
  │   └── ...
  └── ...
```

**Filename format**: `Well{WELL}_Channel{...}_Seq{SEQ}-MaxIP_XY{XY}_{CHANNEL}.tif`
- Example: `WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_405.tif`
- `WELL`: Well identifier (e.g., D10, E10, A01)
- `SEQ`: Sequence number (e.g., 0009, 0010)
- `XY`: XY position number (e.g., 1, 2, 9)
- `CHANNEL`: `405` for Hoechst (nuclear), `488` for MIRO160mer (mitochondrial)
- Files are automatically matched by `(well, sequence, XY)` tuple
- Multiple XY positions and sequences per well are fully supported

## Usage

### Using uv run (recommended)

After running `uv sync`, use `uv run` to execute commands:

```bash
# Basic usage (all dates, default parameters)
uv run edge-spot-pipeline --input inputs/ --output results/

# Process specific date
uv run edge-spot-pipeline --input inputs/ --output results/ --date 20231115

# Use custom parameter configuration
uv run edge-spot-pipeline --input inputs/ --output results/ --config params.json

# Verbose logging
uv run edge-spot-pipeline --input inputs/ --output results/ --verbose
```

### Alternative: Run as module

You can also run the pipeline as a module:

```bash
uv run python -m edge_spot_analyser.pipeline --input inputs/ --output results/
```

### After installation (if using pip)

If you installed with pip, you can run commands directly:

```bash
edge-spot-pipeline --input inputs/ --output results/
```

## Parameter Configuration

You can customize parameters per date using a JSON or YAML configuration file.

### Example `params.json`:

```json
{
  "default": {
    "otsu_correction_factor": 0.45,
    "diameter_min": 30,
    "diameter_max": 100,
    "min_distance": 30,
    "edge_spot_diameter_min": 5,
    "edge_spot_diameter_max": 80,
    "edge_spot_correction_factor": 2.0,
    "inner_expansion": 10,
    "outer_expansion": 15
  },
  "20231115": {
    "otsu_correction_factor": 0.50
  },
  "20231120": {
    "otsu_correction_factor": 0.40,
    "diameter_min": 25
  }
}
```

### Key Parameters

#### Nuclei Segmentation
- **`otsu_correction_factor`** (default: 0.45): Multiplier for Otsu threshold. Lower = more nuclei detected. **This is the main parameter to adjust per date.**
- **`diameter_min`** (default: 30): Minimum nucleus diameter in pixels
- **`diameter_max`** (default: 100): Maximum nucleus diameter in pixels
- **`min_distance`** (default: 30): Minimum distance between nuclei for declumping

#### Edge Spot Detection
- **`edge_spot_correction_factor`** (default: 2.0): Multiplier for Robust Background threshold. Higher = fewer spots detected.
- **`edge_spot_diameter_min`** (default: 5): Minimum spot diameter in pixels
- **`edge_spot_diameter_max`** (default: 80): Maximum spot diameter in pixels

#### Perinuclear Regions
- **`inner_expansion`** (default: 10): Pixels to expand nuclei for inner boundary
- **`outer_expansion`** (default: 15): Pixels to expand for outer boundary (masking)

## Output Files

The pipeline generates CellProfiler-compatible CSV files:

```
results/
  └── 20231115/
      ├── Nuclei.csv              # Nuclear measurements (Hoechst intensity)
      ├── Expand_Nuclei.csv       # Expanded nuclei (MIRO160mer, Gini, moments, etc.)
      ├── Perinuclear_region.csv  # Perinuclear ring measurements
      ├── edge_spots.csv          # Edge spot measurements (if any detected)
      └── Image.csv               # Image-level statistics
```

### Output Measurements

**Nuclei.csv**:
- Basic Hoechst intensity statistics (mean, median, std, integrated, etc.)

**Expand_Nuclei.csv** and **Perinuclear_region.csv**:
- MIRO160mer intensity statistics
- **Gini coefficient** (`GINI_Gini_MIRO160mer`): Inequality measure (0=uniform, 1=concentrated)
- **Mass displacement** (`Intensity_MassDisplacement_MIRO160mer`): Distance between intensity centroid and geometric centroid
- **Moments**: Mean, standard deviation, skewness, kurtosis
- **Edge intensity**: Mean intensity at object border

**edge_spots.csv**:
- MIRO160mer intensity for each detected spot
- Edge intensity

**Image.csv**:
- Global image statistics (total intensity, mean, median, etc.)

## Post-Processing

After generating CSVs, you can run the existing analysis script:

```bash
python ../cellprofiler_output_analyser.py --input results/20231115/
```

This will extract the key metrics:
- **Edge spot fraction**: edge_spot_count / nuclei_count per field of view
- **Mass displacement**: From `Intensity_MassDisplacement_MIRO160mer`
- **Coefficient of Variation (CoV)**: std_intensity / mean_intensity

## Pipeline Architecture

### Module Mapping (from CellProfiler)

| CellProfiler Module | Python Implementation |
|---------------------|----------------------|
| Module 7: IdentifyPrimaryObjects (Nuclei) | `segmentation.segment_nuclei()` |
| Module 8-10: ExpandOrShrinkObjects, MaskImage | `segmentation.create_perinuclear_regions()` |
| Module 11: IdentifyPrimaryObjects (edge_spots) | `segmentation.detect_edge_spots()` |
| Module 13: FilterObjects | `segmentation.filter_edge_spots_by_edge_intensity()` |
| Module 14: IdentifyTertiaryObjects | `segmentation.create_perinuclear_regions()` |
| Modules 12, 15-16, 20: MeasureObjectIntensity | `measurements.measure_all_object_properties()` |
| Module 17: CalculateMoments | `measurements.IntensityMeasurements.measure_moments()` |
| Module 19: CalculateGini | `measurements.IntensityMeasurements.measure_gini()` |

### File Structure

```
edge_spot_analyser/
├── src/
│   └── edge_spot_analyser/      # Python package
│       ├── __init__.py          # Package initialization
│       ├── pipeline.py          # Main pipeline script (CLI entry point)
│       ├── segmentation.py      # Segmentation algorithms (Modules 7-14)
│       ├── measurements.py      # Measurement functions (Modules 12-20)
│       ├── io_utils.py          # File I/O, CSV export
│       └── validate_pipeline.py # Validation script
├── pipeline_files/              # Original CellProfiler pipeline
│   ├── *.cpproj
│   ├── *.json
│   ├── COMPLETE_PIPELINE_REFERENCE.md
│   └── MODULE_11_DEFINITIVE.md
├── inputs/                      # Input TIF files (organized by date)
├── example_CP_output/           # Example CellProfiler outputs for validation
├── pyproject.toml               # Package configuration (uv/pip)
├── params_example.json          # Example parameter configuration
├── QUICKSTART.md                # Quick start guide
└── README.md                    # This file
```

## Algorithm Details

### Nuclei Segmentation (Module 7)
- **Method**: 3-class Otsu thresholding (uses upper threshold)
- **Declumping**: Intensity-based watershed with min_distance=30px
- **Smoothing**: Gaussian (σ=10px) applied to image before thresholding
- **Size filtering**: 30-100 pixel diameter
- **Border removal**: Yes

### Edge Spot Detection (Module 11)
- **Method**: Robust Background thresholding (CellProfiler-specific)
- **Algorithm**:
  1. Remove outliers: lowest 30%, highest 10%
  2. Calculate robust_mean = MEDIAN of remaining pixels
  3. Calculate robust_std = STD of remaining pixels
  4. Threshold = (robust_mean + 6×robust_std) × 2.0
  5. Clip to [0.0035, 1.0]
  6. Apply Gaussian smoothing (σ=1.3) to **binary mask** (not image!)
- **Declumping**: NONE (simple connected component labeling)
- **Size filtering**: 5-80 pixel diameter
- **Border removal**: Yes

### Gini Coefficient
```python
# Exact implementation from CellProfiler
flattened = sort(pixels)
npix = len(pixels)
normalization = abs(mean(pixels)) * npix * (npix - 1)
kernel = (2 * arange(1, npix+1) - npix - 1) * abs(pixels)
gini = sum(kernel) / normalization
```

### Mass Displacement
```python
mass_displacement = distance(intensity_centroid, geometric_centroid)
```

## Troubleshooting

### No nuclei detected
- **Decrease** `otsu_correction_factor` (try 0.35-0.40)
- Check Hoechst image quality and normalization
- Adjust `diameter_min` if nuclei are smaller

### Too many nuclei (false positives)
- **Increase** `otsu_correction_factor` (try 0.50-0.55)
- Increase `min_distance` to separate clumped objects

### No edge spots detected
- **Decrease** `edge_spot_correction_factor` (try 1.5-1.8)
- Check MIRO160mer image quality
- Verify peripheral regions are not over-masked

### Too many edge spots (false positives)
- **Increase** `edge_spot_correction_factor` (try 2.5-3.0)
- Adjust `edge_spot_diameter_min` to exclude small noise

## Validation

To validate against CellProfiler outputs:

```bash
# Run pipeline on example data
uv run edge-spot-pipeline --input example_inputs/ --output validation_output/

# Compare CSVs
uv run edge-spot-validate \
  --cp_output example_CP_output/24111_siRNA/ \
  --py_output validation_output/
```

Or using the module form:

```bash
uv run python -m edge_spot_analyser.pipeline --input example_inputs/ --output validation_output/
uv run python -m edge_spot_analyser.validate_pipeline \
  --cp_output example_CP_output/24111_siRNA/ \
  --py_output validation_output/
```

Expected agreement:
- Object counts: Exact match (±1-2 due to random seed in watershed)
- Intensity measurements: <5% difference
- Spatial metrics (Gini, mass displacement): <5% difference

## References

- Original CellProfiler pipeline: `pipeline_files/231115_combined_pipeline_new_nomito_fixed_moments.cpproj`
- Complete algorithm documentation: `pipeline_files/COMPLETE_PIPELINE_REFERENCE.md`
- Robust Background algorithm: `pipeline_files/MODULE_11_DEFINITIVE.md`
- CellProfiler version: 4.2.4

## License

This code is a port of a CellProfiler pipeline. CellProfiler is licensed under BSD 3-Clause License.

## Authors

- Original CellProfiler pipeline: Christina (parameter fitting and pipeline design)
- Python port: Generated from CellProfiler pipeline v6
