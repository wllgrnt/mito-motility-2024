# Quick Start Guide - Using uv

This guide shows how to get started with the edge-spot-analyser package using `uv`.

## Prerequisites

1. **Install uv** (if not already installed):
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Or on macOS with Homebrew:
```bash
brew install uv
```

## Installation

### Recommended: Use uv sync

```bash
cd edge_spot_analyser

# Sync dependencies and install the package
uv sync

# Or sync with optional dependencies
uv sync --extra yaml     # YAML config support
uv sync --extra dev      # Development tools (ruff)
uv sync --all-extras     # All optional dependencies
```

This creates a virtual environment (`.venv`) and installs all dependencies automatically.

### Alternative: Traditional pip install

If you prefer pip:

```bash
pip install -e .
```

## Running the Pipeline

After running `uv sync`, use `uv run` to execute commands:

### Process all dates

```bash
uv run edge-spot-pipeline --input ../inputs/ --output ../results/
```

### Process a specific date

```bash
uv run edge-spot-pipeline --input ../inputs/ --output ../results/ --date 20231115
```

### Use custom parameters

```bash
uv run edge-spot-pipeline \
  --input ../inputs/ \
  --output ../results/ \
  --config params_example.json \
  --verbose
```

### Validate outputs

```bash
uv run edge-spot-validate \
  --cp_output example_CP_output/24111_siRNA/ \
  --py_output ../results/20231115/
```

## Development Workflow

### Running without installation

You can run as a module without explicit installation:

```bash
uv run python -m edge_spot_analyser.pipeline --input ../inputs/ --output ../results/
```

### Linting

```bash
# Check code
uv run ruff check .

# Auto-fix issues
uv run ruff check --fix .

# Format code
uv run ruff format .
```

### Running tests

```bash
# Run validation on example data
uv run edge-spot-validate \
  --cp_output example_CP_output/24111_siRNA/ \
  --py_output validation_output/
```

## Using as a Python Library

You can also import and use the package in your own Python scripts:

```python
from edge_spot_analyser import (
    Pipeline,
    NucleiSegmentationParams,
    EdgeSpotParams,
    ImageLoader,
    FileDiscovery,
)

# Create pipeline with custom parameters
nuclei_params = NucleiSegmentationParams(
    otsu_correction_factor=0.50,
    diameter_min=30,
    diameter_max=100,
)

pipeline = Pipeline(nuclei_params=nuclei_params)

# Process a single date
pipeline.process_date_batch(
    input_dir="inputs/",
    output_dir="results/",
    date="20231115"
)
```

## Quick Parameter Tuning

If you need to adjust parameters for a specific date:

1. Copy the example config:
```bash
cp params_example.json my_params.json
```

2. Edit `my_params.json`:
```json
{
  "default": {
    "otsu_correction_factor": 0.45
  },
  "20231115": {
    "otsu_correction_factor": 0.50
  }
}
```

3. Run with your config:
```bash
uv run edge-spot-pipeline --input ../inputs/ --output ../results/ --config my_params.json
```

## Troubleshooting

### Command not found: edge-spot-pipeline

Make sure you've run `uv sync` and are using `uv run`:

```bash
# Sync dependencies first
uv sync

# Then run with uv run
uv run edge-spot-pipeline --help
```

### Import errors

If you get import errors, make sure you're in the project directory:

```bash
cd /path/to/mito-motility-2024/edge_spot_analyser
uv run python -m edge_spot_analyser.pipeline --help
```

### Dependencies not installed

If you see missing dependency errors, resync:

```bash
uv sync --reinstall
```

### No nuclei detected

Try lowering the `otsu_correction_factor` in your config file (e.g., from 0.45 to 0.35).

### Too many false positive nuclei

Try increasing the `otsu_correction_factor` (e.g., from 0.45 to 0.55).

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Review the [algorithm details](README.md#algorithm-details)
- Check out the [troubleshooting guide](README.md#troubleshooting)
- Explore the [CellProfiler pipeline reference](pipeline_files/COMPLETE_PIPELINE_REFERENCE.md)
