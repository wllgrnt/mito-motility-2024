"""Tests for core segmentation algorithms."""

import numpy as np
from skimage import morphology

from edge_spot_analyser.segmentation import (
    PerinuclearRegionParams,
    _filter_objects_by_size,
    create_perinuclear_regions,
    detect_edge_spots,
    filter_edge_spots_by_nuclei_proximity,
    robust_background_threshold,
)


class TestRobustBackgroundThreshold:
    def test_known_input(self):
        """Threshold on a known uniform-ish image should be predictable."""
        rng = np.random.default_rng(0)
        # Create image with background ~0.1 and no outliers
        image = rng.normal(loc=0.1, scale=0.01, size=(100, 100)).clip(0, 1).astype(np.float64)

        threshold = robust_background_threshold(
            image,
            lower_outlier_fraction=0.3,
            upper_outlier_fraction=0.1,
            num_deviations=2.0,
            correction_factor=1.0,
            lower_bound=0.0,
            upper_bound=1.0,
        )

        # Median ~0.1, std ~0.01, so threshold ≈ 0.1 + 2*0.01 = 0.12
        assert 0.08 < threshold < 0.20

    def test_all_zeros(self):
        """All-zero image should return lower_bound."""
        image = np.zeros((50, 50))
        threshold = robust_background_threshold(image, lower_bound=0.0035)
        assert threshold == 0.0035

    def test_respects_lower_bound(self):
        """Very dim image threshold should be clipped to lower_bound."""
        image = np.full((50, 50), 0.0001)
        threshold = robust_background_threshold(image, lower_bound=0.0035)
        assert threshold >= 0.0035

    def test_respects_upper_bound(self):
        """Very bright image threshold should be clipped to upper_bound."""
        image = np.ones((50, 50))
        threshold = robust_background_threshold(image, upper_bound=0.5)
        assert threshold <= 0.5


class TestFilterObjectsBySize:
    def test_keeps_correctly_sized(self):
        """Objects within size range are kept."""
        labels = np.zeros((50, 50), dtype=int)
        # Object 1: ~78 pixels (radius 5 disk)
        rr, cc = np.ogrid[:50, :50]
        labels[((rr - 10) ** 2 + (cc - 10) ** 2) <= 25] = 1
        # Object 2: ~12 pixels (radius 2 disk)
        labels[((rr - 30) ** 2 + (cc - 30) ** 2) <= 4] = 2

        # Keep only objects with area 50-200
        result = _filter_objects_by_size(labels, size_min=50, size_max=200)
        assert np.any(result == 1)  # ~78 px, should be kept
        assert not np.any(result == 2)  # ~12 px, should be rejected

    def test_rejects_oversized(self):
        """Objects above size_max are rejected."""
        labels = np.zeros((100, 100), dtype=int)
        labels[10:50, 10:50] = 1  # 1600 pixels
        result = _filter_objects_by_size(labels, size_min=10, size_max=100)
        assert result.max() == 0

    def test_empty_input(self):
        """Empty label image returns empty."""
        labels = np.zeros((20, 20), dtype=int)
        result = _filter_objects_by_size(labels, size_min=10, size_max=100)
        assert result.max() == 0


class TestCreatePerinuclearRegions:
    def test_ring_shape(self):
        """Perinuclear ring should be a donut around each nucleus."""
        # Create a single nucleus at center of 100x100 image
        labels = np.zeros((100, 100), dtype=int)
        rr, cc = np.ogrid[:100, :100]
        labels[((rr - 50) ** 2 + (cc - 50) ** 2) <= 100] = 1  # radius ~10

        params = PerinuclearRegionParams(inner_expansion=10, outer_expansion=15)
        expand_10, expand_15, ring = create_perinuclear_regions(labels, params)

        # Expansion should be larger than original
        assert expand_10.sum() > labels.sum()
        assert expand_15.sum() > expand_10.sum()

        # Ring should not overlap with eroded nucleus
        eroded = morphology.erosion(labels, morphology.disk(1))
        assert not np.any((ring > 0) & (eroded > 0))

        # Ring should be contained within the 10px expansion
        assert np.all((ring > 0) <= (expand_10 > 0))


class TestDetectEdgeSpots:
    def test_blank_image_returns_empty(self):
        """A blank (zero) image should produce no edge spots."""
        blank = np.zeros((100, 100), dtype=np.float64)
        labels = detect_edge_spots(blank)
        assert labels.max() == 0

    def test_detects_bright_spots(self):
        """Bright spots on dark background should be detected."""
        rng = np.random.default_rng(42)
        image = rng.normal(loc=0.05, scale=0.005, size=(200, 200)).clip(0, 1).astype(np.float64)

        # Add bright spots
        image[50:58, 50:58] = 0.8
        image[150:158, 150:158] = 0.8

        labels = detect_edge_spots(image)
        assert labels.max() >= 1


class TestFilterEdgeSpotsByNucleiProximity:
    def test_distant_spots_kept(self):
        """Spots far from nuclei mask should be kept."""
        spots = np.zeros((100, 100), dtype=int)
        spots[80:85, 80:85] = 1  # far from nucleus

        nuclei_mask = np.zeros((100, 100), dtype=int)
        nuclei_mask[10:30, 10:30] = 1  # nucleus region

        result = filter_edge_spots_by_nuclei_proximity(spots, nuclei_mask, expansion_radius=3)
        assert result.max() > 0

    def test_adjacent_spots_removed(self):
        """Spots touching the nuclei mask (after expansion) should be removed."""
        spots = np.zeros((100, 100), dtype=int)
        spots[32:37, 15:20] = 1  # right at edge of nuclei mask

        nuclei_mask = np.zeros((100, 100), dtype=int)
        nuclei_mask[10:32, 10:30] = 1  # nucleus region ends at row 32

        result = filter_edge_spots_by_nuclei_proximity(spots, nuclei_mask, expansion_radius=3)
        assert result.max() == 0

    def test_empty_spots_returns_empty(self):
        """No spots input returns no spots output."""
        spots = np.zeros((50, 50), dtype=int)
        nuclei_mask = np.zeros((50, 50), dtype=int)
        result = filter_edge_spots_by_nuclei_proximity(spots, nuclei_mask, expansion_radius=3)
        assert result.max() == 0
