"""Tests for measurement functions."""

import numpy as np

from edge_spot_analyser.measurements import (
    IntensityMeasurements,
    _calculate_gini_on_pixels,
    measure_edge_spot_burden,
)


class TestGini:
    def test_uniform_array(self):
        """Uniform array should have Gini ≈ 0."""
        pixels = np.full(1000, 0.5)
        gini = _calculate_gini_on_pixels(pixels)
        assert abs(gini) < 0.01

    def test_concentrated(self):
        """All intensity in one pixel should have Gini ≈ 1."""
        pixels = np.zeros(1000)
        pixels[0] = 1.0
        gini = _calculate_gini_on_pixels(pixels)
        assert gini > 0.95

    def test_empty_array(self):
        """Empty array should return 0."""
        gini = _calculate_gini_on_pixels(np.array([]))
        assert gini == 0.0

    def test_all_zeros(self):
        """All-zero array should return 0 (no division error)."""
        gini = _calculate_gini_on_pixels(np.zeros(100))
        assert gini == 0.0

    def test_increasing_inequality(self):
        """More unequal distribution should give higher Gini."""
        uniform = np.full(100, 1.0)
        skewed = np.zeros(100)
        skewed[:10] = 10.0

        gini_uniform = _calculate_gini_on_pixels(uniform)
        gini_skewed = _calculate_gini_on_pixels(skewed)
        assert gini_skewed > gini_uniform


class TestMeasureBasicIntensity:
    def test_known_array(self):
        """Verify mean, sum, etc. on a known labeled image."""
        labels = np.zeros((20, 20), dtype=int)
        labels[5:10, 5:10] = 1  # 25 pixels

        intensity = np.zeros((20, 20), dtype=np.float64)
        intensity[5:10, 5:10] = 0.4

        df = IntensityMeasurements.measure_basic_intensity(labels, intensity)

        assert len(df) == 1
        assert abs(df["Intensity_MeanIntensity"].iloc[0] - 0.4) < 1e-6
        assert abs(df["Intensity_IntegratedIntensity"].iloc[0] - 10.0) < 1e-6
        assert abs(df["Intensity_MinIntensity"].iloc[0] - 0.4) < 1e-6
        assert abs(df["Intensity_MaxIntensity"].iloc[0] - 0.4) < 1e-6

    def test_empty_labels(self):
        """No objects should return empty dataframe."""
        labels = np.zeros((20, 20), dtype=int)
        intensity = np.ones((20, 20), dtype=np.float64)
        df = IntensityMeasurements.measure_basic_intensity(labels, intensity)
        assert len(df) == 0

    def test_multiple_objects(self):
        """Multiple objects should each get their own row."""
        labels = np.zeros((30, 30), dtype=int)
        labels[2:7, 2:7] = 1
        labels[15:20, 15:20] = 2

        intensity = np.zeros((30, 30), dtype=np.float64)
        intensity[2:7, 2:7] = 0.3
        intensity[15:20, 15:20] = 0.6

        df = IntensityMeasurements.measure_basic_intensity(labels, intensity)
        assert len(df) == 2

        obj1 = df[df["ObjectNumber"] == 1].iloc[0]
        obj2 = df[df["ObjectNumber"] == 2].iloc[0]
        assert abs(obj1["Intensity_MeanIntensity"] - 0.3) < 1e-6
        assert abs(obj2["Intensity_MeanIntensity"] - 0.6) < 1e-6


class TestMeasureEdgeSpotBurden:
    def test_known_geometry(self):
        """Known geometry should produce expected ratios."""
        # Create a simple scenario
        edge_spots = np.zeros((100, 100), dtype=int)
        edge_spots[80:90, 80:90] = 1  # 100 pixels of edge spots

        perinuclear = np.zeros((100, 100), dtype=int)
        perinuclear[10:30, 10:30] = 1  # 400 pixels of perinuclear

        # Uniform intensity image
        miro = np.full((100, 100), 0.5, dtype=np.float64)

        metrics = measure_edge_spot_burden(
            edge_spots_labels=edge_spots,
            miro_image=miro,
            perinuclear_region_labels=perinuclear,
            nuclei_count=5,
        )

        # Edge spot intensity: 100 pixels * 0.5 = 50.0
        assert abs(metrics["edge_spot_intensity_total"] - 50.0) < 1e-6

        # Per nucleus: 50.0 / 5 = 10.0
        assert abs(metrics["edge_spot_intensity_per_nucleus"] - 10.0) < 1e-6

        # Fraction of total: 50.0 / (100*100*0.5) = 50/5000 = 0.01
        assert abs(metrics["edge_fraction_of_total_miro"] - 0.01) < 1e-6

        # Edge/perinuclear ratio: 50.0 / (400*0.5) = 50/200 = 0.25
        assert abs(metrics["edge_to_perinuclear_ratio"] - 0.25) < 1e-6

    def test_no_edge_spots(self):
        """No edge spots should give zero metrics."""
        edge_spots = np.zeros((50, 50), dtype=int)
        perinuclear = np.zeros((50, 50), dtype=int)
        perinuclear[10:20, 10:20] = 1
        miro = np.full((50, 50), 0.5, dtype=np.float64)

        metrics = measure_edge_spot_burden(edge_spots, miro, perinuclear, nuclei_count=3)
        assert metrics["edge_spot_intensity_total"] == 0.0
        assert metrics["edge_spot_intensity_per_nucleus"] == 0.0

    def test_zero_nuclei(self):
        """Zero nuclei count should return zero per-nucleus metric."""
        edge_spots = np.zeros((50, 50), dtype=int)
        edge_spots[5:10, 5:10] = 1
        perinuclear = np.zeros((50, 50), dtype=int)
        miro = np.full((50, 50), 0.5, dtype=np.float64)

        metrics = measure_edge_spot_burden(edge_spots, miro, perinuclear, nuclei_count=0)
        assert metrics["edge_spot_intensity_per_nucleus"] == 0.0
