# mito_motility_2024

[![DOI](https://zenodo.org/badge/836738016.svg)](https://zenodo.org/doi/10.5281/zenodo.13150979)

Scripts used in Christina Gladkova's 2024 paper on mitochondrial motility.

Contents:
- Gini plugin: used to generate Gini coefficients with CellProfiler, as part of a pipeline including e.g. image segmentation.
- cellprofiler_output_analyser.py: Given a fixed input folder structure with CellProfiler CSV files, extracts the relevant data and stacks
    into a tidier columnar format.
- trackmate_analyser.ipynb: rotates all single-particle tracks to be a consistent direction and extracts e.g. the distribution of speeds.
