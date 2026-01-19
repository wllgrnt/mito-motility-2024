# Edge Spot Analyser

Python port of a CellProfiler pipeline for analyzing mitochondrial distribution around nuclei.

## Overview

This pipeline analyzes the distribution of mitochondria (MIRO160mer) around nuclei in microscopy images. It identifies bright peripheral mitochondrial spots and quantifies their spatial distribution using metrics like Gini coefficient, mass displacement, and coefficient of variation.

### What it does:

1. **Segments nuclei** from Hoechst (405nm) channel using 3-class Otsu thresholding
2. **Defines perinuclear regions** (ring between eroded nuclei and 10px expansion)
3. **Detects bright peripheral mitochondrial spots** in MIRO160mer (488nm) channel
4. **Quantifies distribution metrics**: Gini coefficient, coefficient of variation (CoV), mass displacement, moments (skewness, kurtosis)
5. **Aggregates results** per date with edge spot fractions and summary statistics
6. **Generates figure tables** combining data across dates for publication figures

## Quick Start

### Prerequisites

Install uv if not already installed:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# Or on macOS: brew install uv
```

### Installation

```bash
cd edge_spot_analyser
uv sync
```

### Run the Pipeline

```bash
# Process all dates (includes per-date aggregation)
uv run edge-spot-pipeline -i inputs/ -o results/

# Process specific date with parallel workers
uv run edge-spot-pipeline -i inputs/ -o results/ -d 241223 -w 10

# Use custom parameters from Excel config
uv run edge-spot-pipeline -i inputs/ -o results/ -c "christina obligations.xlsx"

# Skip aggregation (only run image processing)
uv run edge-spot-pipeline -i inputs/ -o results/ --skip-aggregate
```

### Run Aggregation Separately

```bash
# Per-date aggregation (edge spot fractions, Gini, etc.)
uv run edge-spot-aggregate -i results/

# Figure aggregation (combine across dates for publication figures)
uv run edge-spot-figure-aggregate -r results/ -c "christina obligations.xlsx"
```

## Input Data Format

The pipeline expects TIF files organized in date-based subdirectories:

```
inputs/
  ├── 241223/
  │   ├── WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_405.tif
  │   ├── WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_488.tif
  │   └── ...
  └── 241217/
      └── ...
```

**Filename format**: `Well{WELL}_Channel{...}_Seq{SEQ}[-MaxIP]_XY{XY}_{CHANNEL}.tif`
- `WELL`: Well identifier (e.g., B02, D10)
- `SEQ`: Sequence number (e.g., 0000, 0009)
- `-MaxIP`: Optional (maximum intensity projection indicator)
- `XY`: XY position number (1-9)
- `CHANNEL`: `405` for Hoechst (nuclear), `488` for MIRO160mer (mitochondrial)

Files are automatically matched by `(well, sequence, XY)` tuple.

## Parameter Configuration

Parameters can be configured via JSON file or Excel spreadsheet.

### JSON Configuration

```json
{
  "default": {
    "otsu_correction_factor": 0.45,
    "diameter_min": 30,
    "diameter_max": 100,
    "min_distance": 30,
    "edge_spot_diameter_min": 5,
    "edge_spot_diameter_max": 80,
    "edge_spot_correction_factor": 2.0
  },
  "241223": {
    "otsu_correction_factor": 0.50
  }
}
```

### Excel Configuration

The pipeline can read parameters from an Excel file with an `Otsu_params` sheet containing per-date settings and exclusions.

### Key Parameters

**Nuclei Segmentation**
- **`otsu_correction_factor`** (default: 0.45): Multiplier for Otsu threshold. Lower = more nuclei detected.
- **`diameter_min`** / **`diameter_max`** (default: 30/100): Nucleus diameter range in pixels
- **`min_distance`** (default: 30): Minimum distance between nuclei for declumping

**Edge Spot Detection**
- **`edge_spot_correction_factor`** (default: 2.0): Multiplier for threshold. Higher = fewer spots detected.
- **`edge_spot_diameter_min`** / **`edge_spot_diameter_max`** (default: 5/80): Spot diameter range in pixels

## Output Files

### Per-Date Output (from pipeline)

```
results/
  └── 241223/
      ├── Nuclei.csv                    # Nuclear measurements (Hoechst intensity)
      ├── Expand_Nuclei.csv             # Expanded nuclei (MIRO160mer, Gini, moments)
      ├── Perinuclear_region.csv        # Perinuclear ring measurements
      ├── edge_spots.csv                # Edge spot measurements
      ├── Image.csv                     # Image-level statistics
      ├── All_measurements.csv          # Combined summary
      ├── edge_spot_fraction_static.csv # Edge spot fractions per well
      ├── GINI_Gini_MIRO160mer_fov_median_static.csv  # Gini per FOV
      └── ...
```

### Figure Output (from figure aggregation)

```
results/
  └── figures/
      ├── Fig1D_edge_spot_raw.csv        # Raw values for Figure 1D
      ├── Fig1D_edge_spot_normalized.csv # Normalized (per-date, control=1.0)
      ├── Fig1D_gini_raw.csv
      ├── Fig1D_gini_normalized.csv
      └── ...
```

### Key Measurements

- **Gini coefficient** (`GINI_Gini_MIRO160mer`): Inequality measure (0=uniform, 1=concentrated)
- **Mass displacement** (`Intensity_MassDisplacement_MIRO160mer`): Distance between intensity and geometric centroids
- **Moments**: Mean, standard deviation, skewness, kurtosis
- **Edge spot fraction**: Proportion of nuclei with associated edge spots

## Validation

Compare Python output with CellProfiler reference:

```bash
uv run python -m edge_spot_analyser.validate_pipeline \
  --cp_output example_CP_output/241223_siRNA_motor_JNK_wt_DRH/ \
  --py_output results/241223/
```

## Algorithm Details

### Nuclei Segmentation (Module 7)
- **Method**: 3-class Otsu thresholding (uses upper threshold)
- **Smoothing**: Gaussian (σ=1.0, from threshold smoothing scale=2.0)
- **Declumping**: Intensity-based watershed with min_distance=30px
- **Size filtering**: 30-100 pixel diameter
- **Border removal**: Yes

### Perinuclear Region (Modules 8, 9, 14)
- **Expand_Nuclei**: Nuclei expanded by 10 pixels (for measurements)
- **Expand_Nuclei_for_mask**: Nuclei expanded by 15 pixels (for masking only)
- **Perinuclear_region**: Ring between eroded nuclei (1px erosion) and 10px expansion

### Edge Spot Detection (Module 11 - Robust Background)
1. Remove outliers: lowest 30%, highest 10%
2. Calculate robust_mean = MEDIAN of remaining pixels
3. Calculate robust_std = STD of remaining pixels
4. Threshold = (robust_mean + 6×robust_std) × correction_factor
5. Clip to [0.0035, 1.0]
6. Apply Gaussian smoothing (σ=1.3) to binary mask
7. Simple connected component labeling (no declumping)

### Gini Coefficient (Module 19)
```python
flattened = sort(pixels)
npix = len(pixels)
normalization = abs(mean(pixels)) * npix * (npix - 1)
kernel = (2 * arange(1, npix+1) - npix - 1) * abs(pixels)
gini = sum(kernel) / normalization
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No nuclei detected | Decrease `otsu_correction_factor` (try 0.35-0.40) |
| Too many nuclei | Increase `otsu_correction_factor` (try 0.50-0.55) |
| No edge spots | Decrease `edge_spot_correction_factor` (try 1.5-1.8) |
| Too many edge spots | Increase `edge_spot_correction_factor` (try 2.5-3.0) |
| Images not found | Check filename format matches expected pattern |

## File Structure

```
edge_spot_analyser/
├── src/edge_spot_analyser/
│   ├── pipeline.py            # Main pipeline (CLI: edge-spot-pipeline)
│   ├── segmentation.py        # Segmentation algorithms
│   ├── measurements.py        # Measurement functions
│   ├── io_utils.py            # File I/O, CSV export
│   ├── aggregation.py         # Per-date aggregation (CLI: edge-spot-aggregate)
│   ├── figure_aggregation.py  # Cross-date figure tables (CLI: edge-spot-figure-aggregate)
│   └── validate_pipeline.py   # Validation against CellProfiler
├── pipeline_files/            # Original CellProfiler pipeline JSON
├── inputs/                    # Input images (by date)
├── results/                   # Output CSVs (by date)
├── pyproject.toml             # Package configuration
└── README.md
```

## License

MIT License. See [LICENSE](../LICENSE) for details.
