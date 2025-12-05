# Edge Spot Analyser

Python port of a CellProfiler pipeline for analyzing mitochondrial distribution around nuclei.

## Overview

This pipeline analyzes the distribution of mitochondria (MIRO160mer) around nuclei in microscopy images. It identifies bright peripheral mitochondrial spots and quantifies their spatial distribution using metrics like Gini coefficient, mass displacement, and coefficient of variation.

### What it does:

1. **Segments nuclei** from Hoechst (405nm) channel using 3-class Otsu thresholding
2. **Defines perinuclear regions** (10-15 pixel rings around nuclei)
3. **Detects bright peripheral mitochondrial spots** in MIRO160mer (488nm) channel
4. **Quantifies distribution metrics**: Gini coefficient, coefficient of variation (CoV), mass displacement, moments (skewness, kurtosis)
5. **Exports CellProfiler-compatible CSV files** for downstream analysis

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
# Process all dates
uv run python -m edge_spot_analyser.pipeline --input inputs/ --output results/

# Process specific date with parallel workers
uv run python -m edge_spot_analyser.pipeline --input inputs/ --output results/ --date 241223 --workers 10

# Use custom parameters
uv run python -m edge_spot_analyser.pipeline --input inputs/ --output results/ --config params.json
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

**Filename format**: `Well{WELL}_Channel{...}_Seq{SEQ}-MaxIP_XY{XY}_{CHANNEL}.tif`
- `WELL`: Well identifier (e.g., B02, D10)
- `SEQ`: Sequence number (e.g., 0000, 0009)
- `XY`: XY position number (1-9)
- `CHANNEL`: `405` for Hoechst (nuclear), `488` for MIRO160mer (mitochondrial)

Files are automatically matched by `(well, sequence, XY)` tuple.

## Parameter Configuration

Customize parameters per date using a JSON configuration file:

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
  "241223": {
    "otsu_correction_factor": 0.50
  }
}
```

### Key Parameters

**Nuclei Segmentation**
- **`otsu_correction_factor`** (default: 0.45): Multiplier for Otsu threshold. Lower = more nuclei detected.
- **`diameter_min`** / **`diameter_max`** (default: 30/100): Nucleus diameter range in pixels
- **`min_distance`** (default: 30): Minimum distance between nuclei for declumping

**Edge Spot Detection**
- **`edge_spot_correction_factor`** (default: 2.0): Multiplier for threshold. Higher = fewer spots detected.
- **`edge_spot_diameter_min`** / **`edge_spot_diameter_max`** (default: 5/80): Spot diameter range in pixels

**Perinuclear Regions**
- **`inner_expansion`** (default: 10): Pixels to expand nuclei for inner boundary
- **`outer_expansion`** (default: 15): Pixels to expand for outer boundary

## Output Files

```
results/
  └── 241223/
      ├── Nuclei.csv              # Nuclear measurements (Hoechst intensity)
      ├── Expand_Nuclei.csv       # Expanded nuclei (MIRO160mer, Gini, moments)
      ├── Perinuclear_region.csv  # Perinuclear ring measurements
      ├── edge_spots.csv          # Edge spot measurements
      ├── Image.csv               # Image-level statistics
      └── All_measurements.csv    # Combined summary
```

### Key Measurements

- **Gini coefficient** (`GINI_Gini_MIRO160mer`): Inequality measure (0=uniform, 1=concentrated)
- **Mass displacement** (`Intensity_MassDisplacement_MIRO160mer`): Distance between intensity and geometric centroids
- **Moments**: Mean, standard deviation, skewness, kurtosis
- **Edge intensity**: Mean intensity at object border

## Validation

Compare Python output with CellProfiler reference:

```bash
uv run python -m edge_spot_analyser.validate_pipeline \
  --cp_output example_CP_output/241223_siRNA_motor_JNK_wt_DRH/ \
  --py_output results/241223/
```

## Algorithm Details

### Nuclei Segmentation
- **Method**: 3-class Otsu thresholding (uses upper threshold)
- **Smoothing**: Gaussian (σ=1.0, from threshold smoothing scale=2.0)
- **Declumping**: Intensity-based watershed with min_distance=30px
- **Size filtering**: 30-100 pixel diameter
- **Border removal**: Yes

### Edge Spot Detection (Robust Background)
1. Remove outliers: lowest 30%, highest 10%
2. Calculate robust_mean = MEDIAN of remaining pixels
3. Calculate robust_std = STD of remaining pixels
4. Threshold = (robust_mean + 6×robust_std) × correction_factor
5. Clip to [0.0035, 1.0]
6. Apply Gaussian smoothing (σ=1.3) to binary mask
7. Simple connected component labeling (no declumping)

### Gini Coefficient
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

## File Structure

```
edge_spot_analyser/
├── src/edge_spot_analyser/
│   ├── pipeline.py          # Main pipeline (CLI entry point)
│   ├── segmentation.py      # Segmentation algorithms
│   ├── measurements.py      # Measurement functions
│   ├── io_utils.py          # File I/O, CSV export
│   └── validate_pipeline.py # Validation against CellProfiler
├── pipeline_files/          # Original CellProfiler pipeline
├── example_CP_output/       # Reference outputs for validation
├── inputs/                  # Input images (by date)
├── results/                 # Output CSVs (by date)
├── pyproject.toml           # Package configuration
└── params_example.json      # Example parameter config
```

## License

MIT License. See [LICENSE](../LICENSE) for details.
