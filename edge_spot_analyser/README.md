# Edge Spot Analyser

A re-implementation of the CellProfiler pipeline initially used to generate the figures in https://www.biorxiv.org/content/10.1101/2024.09.13.612963v1.

## Citation

If you use this software, please cite:

> [Paper citation placeholder - to be updated upon publication]

## Raw Data Availability

The raw microscopy data (~240GB of TIF files) is available on Zenodo:

> **DOI:** [Zenodo DOI placeholder - to be filled after data upload]

After downloading, extract the data into the `inputs/` directory to match the expected structure described below.

The prior pipeline's spec can be found in `pipeline_files/` - this output was then analysed using `cellprofiler_output_analyser.py` in the repo root. During the revision process we needed to re-run the CellProfiler pipeline with different parameters, so it made since to port our pipeline into Python, and consolidate all the analysis work into one location.


## Overview

The original pipeline's steps were:
- Load a set of tif files for a set of dates (including multiple well numbers per date, and multiple xy positions per well, for each channel).
- Assign files ending in 405.tif (i.e 405nm) to the Hoechst channel (nuclear staining), and files ending in 488.tif to the MIRO160mer channel (a synthetic cargo used as a model for mitochondrial location, see paper for details).
- Use 3-class Otsu thresholding to segment nuclei in the Hoechst channel.
- Expand the nuclear mask by 10 pixels to define a perinuclear region.
- Expand the nuclear mask by 15 pixels to define a mask for excluding the perinuclear region.
- Identify 'edge spots' - bright points in the image outside the perinuclear region, using Robust Background thresholding. Measure the edge spot intensities. Count the number of such objects.
- Measure the perinuclear region intensities (in the 10 pixel donut around each nucleus). Compute moments on the perinuclear region, calculate the Gini coefficient of the perinuclear intensities.
- Export to spreadsheet.
- (postprocessing, in `cellprofiler_output_analyser.py`): Generate the median Gini and `edge spot count/ nuclei count` per (Date, WellNumber, XY) and stacks them into columns for each wellnumber. From there we can generate the paper's figures for each set of conditions.


This repo replicates the above steps, with the following differences:

- A smoothing module before the Otsu step.
- Automatically generates the necessary figures from a spreadsheet mapping condition to (date,wellnumber).


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
uv run edge-spot-pipeline -i inputs/ -o results/ -c config.xlsx

# Skip aggregation (only run image processing)
uv run edge-spot-pipeline -i inputs/ -o results/ --skip-aggregate
```

### Run Aggregation Separately

```bash
# Per-date aggregation (edge spot fractions, Gini, etc.)
uv run edge-spot-aggregate -i results/

# Figure aggregation (combine across dates for publication figures)
uv run edge-spot-figure-aggregate -r results/ -c config.xlsx
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
- **`edge_spot_diameter_min`** / **`edge_spot_diameter_max`** (default: 3/80): Spot diameter range in pixels

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
4. Threshold = (robust_mean + 2×robust_std) × correction_factor
5. Clip to [0.0035, 1.0]
6. Apply Gaussian smoothing (σ=0.65) to image before thresholding
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
│   └── plotting.py            # Figure plotting (CLI: edge-spot-plot)
├── pipeline_files/            # Original CellProfiler pipeline JSON
├── inputs/                    # Input images (by date)
├── results/                   # Output CSVs (by date)
├── pyproject.toml             # Package configuration
└── README.md
```

## License

MIT License. See [LICENSE](LICENSE) for details.