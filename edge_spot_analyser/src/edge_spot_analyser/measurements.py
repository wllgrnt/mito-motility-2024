"""
Measurements module for the CellProfiler pipeline port.

This module implements measurement functions from the CellProfiler pipeline:
- Module 12: MeasureObjectIntensity (edge_spots)
- Module 15: MeasureObjectIntensity (Nuclei - Hoechst)
- Module 16: MeasureObjectIntensity (Expand_Nuclei, Perinuclear_region - MIRO160mer)
- Module 17: CalculateMoments (distribution statistics)
- Module 19: CalculateGini (Gini coefficient)
- Module 20: MeasureImageIntensity (global statistics)

Author: Generated from CellProfiler pipeline 231115_combined_pipeline_new_nomito_fixed_moments
Pipeline Version: v6
"""

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from skimage import measure, morphology


class IntensityMeasurements:
    """
    Calculate intensity-based measurements for labeled objects.

    This class provides all intensity measurements from CellProfiler modules,
    including basic statistics, edge intensity, and spatial distribution metrics.
    """

    @staticmethod
    def measure_basic_intensity(
        labels: np.ndarray, intensity_image: np.ndarray, prefix: str = "Intensity"
    ) -> pd.DataFrame:
        """
        Measure basic intensity statistics for each labeled object.

        Implements intensity measurements from MeasureObjectIntensity modules.

        Measurements include:
        - IntegratedIntensity: Sum of all pixel intensities
        - MeanIntensity: Average intensity
        - MedianIntensity: Median intensity
        - StdIntensity: Standard deviation
        - MinIntensity: Minimum intensity
        - MaxIntensity: Maximum intensity
        - LowerQuartileIntensity: 25th percentile
        - UpperQuartileIntensity: 75th percentile
        - MADIntensity: Median Absolute Deviation

        Args:
            labels: Labeled image where each object has a unique integer label
            intensity_image: 2D float array of intensity values
            prefix: Column name prefix (e.g., "Intensity")

        Returns:
            DataFrame with one row per object, columns for each measurement
        """
        if labels.max() == 0:
            # No objects, return empty dataframe
            return pd.DataFrame()

        # Use scikit-image regionprops for basic measurements
        props = measure.regionprops_table(
            labels, intensity_image=intensity_image, properties=["label", "area"]
        )

        df = pd.DataFrame(props)
        df.rename(columns={"label": "ObjectNumber"}, inplace=True)

        # Calculate intensity statistics for each object
        measurements = []

        for region in measure.regionprops(labels, intensity_image=intensity_image):
            pixels = region.intensity_image[region.image]  # Only pixels in this object

            meas = {
                "ObjectNumber": region.label,
                f"{prefix}_IntegratedIntensity": pixels.sum(),
                f"{prefix}_MeanIntensity": pixels.mean(),
                f"{prefix}_MedianIntensity": np.median(pixels),
                f"{prefix}_StdIntensity": pixels.std(),
                f"{prefix}_MinIntensity": pixels.min(),
                f"{prefix}_MaxIntensity": pixels.max(),
                f"{prefix}_LowerQuartileIntensity": np.percentile(pixels, 25),
                f"{prefix}_UpperQuartileIntensity": np.percentile(pixels, 75),
                f"{prefix}_MADIntensity": np.median(np.abs(pixels - np.median(pixels))),
            }
            measurements.append(meas)

        measurements_df = pd.DataFrame(measurements)

        # Merge with basic properties
        df = df.merge(measurements_df, on="ObjectNumber", how="left")

        return df

    @staticmethod
    def measure_edge_intensity(
        labels: np.ndarray, intensity_image: np.ndarray, dilation_radius: int = 1
    ) -> pd.Series:
        """
        Measure mean intensity at the edge of each object.

        Edge is defined as the border region created by dilating the object
        and subtracting the original object.

        Args:
            labels: Labeled image
            intensity_image: 2D float array of intensity values
            dilation_radius: Number of pixels to dilate (default 1)

        Returns:
            Series mapping ObjectNumber to MeanIntensityEdge
        """
        edge_intensities = {}

        for region in measure.regionprops(labels):
            # Get object mask
            mask = labels == region.label

            # Dilate by dilation_radius
            dilated_mask = morphology.binary_dilation(mask, morphology.disk(dilation_radius))

            # Edge is dilated - original
            edge_mask = dilated_mask & ~mask

            # Calculate mean edge intensity
            if edge_mask.sum() > 0:
                edge_intensity = intensity_image[edge_mask].mean()
            else:
                edge_intensity = 0.0

            edge_intensities[region.label] = edge_intensity

        return pd.Series(edge_intensities, name="Intensity_MeanIntensityEdge")

    @staticmethod
    def measure_mass_displacement(labels: np.ndarray, intensity_image: np.ndarray) -> pd.Series:
        """
        Measure mass displacement for each object.

        Mass displacement is the distance between the geometric centroid
        and the intensity-weighted centroid (center of mass).

        A high mass displacement indicates asymmetric intensity distribution.

        Args:
            labels: Labeled image
            intensity_image: 2D float array of intensity values

        Returns:
            Series mapping ObjectNumber to MassDisplacement
        """
        mass_displacements = {}

        for region in measure.regionprops(labels, intensity_image=intensity_image):
            # Geometric centroid (y, x)
            centroid = np.array(region.centroid)

            # Intensity-weighted centroid (center of mass)
            weighted_centroid = np.array(region.centroid_weighted)

            # Euclidean distance
            mass_displacement = np.linalg.norm(centroid - weighted_centroid)

            mass_displacements[region.label] = mass_displacement

        return pd.Series(mass_displacements, name="Intensity_MassDisplacement")

    @staticmethod
    def measure_gini(labels: np.ndarray, intensity_image: np.ndarray) -> pd.Series:
        """
        Calculate Gini coefficient for each object.

        The Gini coefficient measures inequality in intensity distribution:
        - 0 = perfectly uniform (all pixels have same intensity)
        - 1 = perfectly unequal (all intensity in one pixel)

        High Gini indicates clustered/punctate distribution.

        Algorithm from CellProfiler calculate_gini.py:
            flattened = sort(pixels)
            npix = len(pixels)
            normalization = abs(mean(pixels)) * npix * (npix - 1)
            kernel = (2 * arange(1, npix+1) - npix - 1) * abs(pixels)
            gini = sum(kernel) / normalization

        Args:
            labels: Labeled image
            intensity_image: 2D float array of intensity values

        Returns:
            Series mapping ObjectNumber to Gini coefficient
        """
        gini_coefficients = {}

        for region in measure.regionprops(labels, intensity_image=intensity_image):
            pixels = region.intensity_image[region.image]
            gini = _calculate_gini_on_pixels(pixels)
            gini_coefficients[region.label] = gini

        return pd.Series(gini_coefficients, name="GINI_Gini")

    @staticmethod
    def measure_moments(labels: np.ndarray, intensity_image: np.ndarray) -> pd.DataFrame:
        """
        Calculate distribution moments for each object.

        Implements CellProfiler Module 17 (CalculateMoments).

        Measurements:
        - Mean: Average intensity (1st moment)
        - Standard Deviation: Spread (2nd moment)
        - Skewness: Asymmetry (3rd moment, Fisher's definition)
        - Kurtosis: Tailedness (4th moment, Fisher's definition)

        Args:
            labels: Labeled image
            intensity_image: 2D float array of intensity values

        Returns:
            DataFrame with columns for Mean, StdDev, Skewness, Kurtosis
        """
        moments_list = []

        for region in measure.regionprops(labels, intensity_image=intensity_image):
            pixels = region.intensity_image[region.image]

            # Calculate moments using scipy.stats
            # Fisher=True means excess kurtosis (normal distribution = 0)
            moments = {
                "ObjectNumber": region.label,
                "Moments_Mean": pixels.mean(),
                "Moments_StandardDeviation": pixels.std(),
                "Moments_Skewness": stats.skew(pixels),
                "Moments_Kurtosis": stats.kurtosis(pixels, fisher=True),
            }
            moments_list.append(moments)

        return pd.DataFrame(moments_list)


def _calculate_gini_on_pixels(pixels: np.ndarray) -> float:
    """
    Calculate Gini coefficient on a 1D array of pixel intensities.

    This is the exact implementation from CellProfiler's calculate_gini.py plugin.

    Args:
        pixels: 1D array of pixel intensities

    Returns:
        Gini coefficient (0-1)
    """
    flattened = np.sort(np.ravel(pixels))
    npix = np.size(flattened)

    if npix == 0:
        return 0.0

    mean_intensity = np.abs(np.mean(flattened))

    if mean_intensity == 0:
        return 0.0

    normalization = mean_intensity * npix * (npix - 1)

    if normalization == 0:
        return 0.0

    kernel = (2.0 * np.arange(1, npix + 1) - npix - 1) * np.abs(flattened)
    gini = np.sum(kernel) / normalization

    return gini


class ImageMeasurements:
    """Calculate image-level (global) intensity statistics."""

    @staticmethod
    def measure_image_intensity(
        intensity_image: np.ndarray, mask: np.ndarray | None = None
    ) -> dict[str, float]:
        """
        Measure global intensity statistics for an entire image.

        Implements CellProfiler Module 20 (MeasureImageIntensity).

        Args:
            intensity_image: 2D float array of intensity values
            mask: Optional binary mask (only measure within mask)

        Returns:
            Dictionary of image-level measurements
        """
        if mask is not None:
            pixels = intensity_image[mask]
        else:
            pixels = intensity_image.flatten()

        # Remove zeros (background)
        pixels = pixels[pixels > 0]

        if len(pixels) == 0:
            return {
                "Image_Intensity_TotalIntensity": 0.0,
                "Image_Intensity_MeanIntensity": 0.0,
                "Image_Intensity_MedianIntensity": 0.0,
                "Image_Intensity_StdIntensity": 0.0,
                "Image_Intensity_MinIntensity": 0.0,
                "Image_Intensity_MaxIntensity": 0.0,
            }

        return {
            "Image_Intensity_TotalIntensity": pixels.sum(),
            "Image_Intensity_MeanIntensity": pixels.mean(),
            "Image_Intensity_MedianIntensity": np.median(pixels),
            "Image_Intensity_StdIntensity": pixels.std(),
            "Image_Intensity_MinIntensity": pixels.min(),
            "Image_Intensity_MaxIntensity": pixels.max(),
        }


def measure_all_object_properties(
    labels: np.ndarray,
    intensity_image: np.ndarray,
    channel_name: str,
    include_edge_intensity: bool = False,
    include_mass_displacement: bool = False,
    include_gini: bool = False,
    include_moments: bool = False,
) -> pd.DataFrame:
    """
    Measure all properties for a set of labeled objects.

    This is a convenience function that combines all measurements into a single DataFrame.

    Args:
        labels: Labeled image
        intensity_image: 2D float array of intensity values
        channel_name: Name of the channel (e.g., "MIRO160mer", "Hoechst")
        include_edge_intensity: Whether to calculate edge intensity
        include_mass_displacement: Whether to calculate mass displacement
        include_gini: Whether to calculate Gini coefficient
        include_moments: Whether to calculate moments

    Returns:
        DataFrame with all measurements, one row per object
    """
    if labels.max() == 0:
        return pd.DataFrame()

    # Start with basic intensity measurements
    prefix = "Intensity"
    df = IntensityMeasurements.measure_basic_intensity(labels, intensity_image, prefix=prefix)

    # Add optional measurements
    if include_edge_intensity:
        edge_intensity = IntensityMeasurements.measure_edge_intensity(labels, intensity_image)
        df[f"{prefix}_MeanIntensityEdge_{channel_name}"] = df["ObjectNumber"].map(edge_intensity)

    if include_mass_displacement:
        mass_disp = IntensityMeasurements.measure_mass_displacement(labels, intensity_image)
        df[f"{prefix}_MassDisplacement_{channel_name}"] = df["ObjectNumber"].map(mass_disp)

    if include_gini:
        gini = IntensityMeasurements.measure_gini(labels, intensity_image)
        df[f"GINI_Gini_{channel_name}"] = df["ObjectNumber"].map(gini)

    if include_moments:
        moments = IntensityMeasurements.measure_moments(labels, intensity_image)
        moments_renamed = moments.rename(
            columns={
                "Moments_Mean": f"Moments_Mean_{channel_name}",
                "Moments_StandardDeviation": f"Moments_StandardDeviation_{channel_name}",
                "Moments_Skewness": f"Moments_Skewness_{channel_name}",
                "Moments_Kurtosis": f"Moments_Kurtosis_{channel_name}",
            }
        )
        df = df.merge(moments_renamed, on="ObjectNumber", how="left")

    # Rename basic intensity columns to include channel name
    basic_cols = [
        "IntegratedIntensity",
        "MeanIntensity",
        "MedianIntensity",
        "StdIntensity",
        "MinIntensity",
        "MaxIntensity",
        "LowerQuartileIntensity",
        "UpperQuartileIntensity",
        "MADIntensity",
    ]

    rename_dict = {f"{prefix}_{col}": f"{prefix}_{col}_{channel_name}" for col in basic_cols}
    df.rename(columns=rename_dict, inplace=True)

    return df


def combine_measurements_for_export(
    nuclei_labels: np.ndarray,
    expand_nuclei_labels: np.ndarray,
    perinuclear_region_labels: np.ndarray,
    edge_spots_labels: np.ndarray,
    hoechst_image: np.ndarray,
    miro_image: np.ndarray,
    metadata: dict[str, Any],
    masked_miro_image: np.ndarray | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Combine all measurements into CellProfiler-compatible DataFrames.

    This function orchestrates all measurements and organizes them into the
    same structure as CellProfiler's CSV outputs.

    Args:
        nuclei_labels: Labeled nuclei
        expand_nuclei_labels: Expanded nuclei (10px)
        perinuclear_region_labels: Perinuclear ring region
        edge_spots_labels: Detected edge spots
        hoechst_image: Hoechst channel image
        miro_image: MIRO160mer channel image (for Expand_Nuclei/Perinuclear measurements)
        metadata: Dictionary with 'ImageNumber', 'FileName_Hoechst', etc.
        masked_miro_image: Masked MIRO160mer image for edge_spots measurements
            (Subtract_perinucleus_Miro160mer in CellProfiler). If None, uses miro_image.

    Returns:
        Dictionary mapping output names to DataFrames:
        - 'Nuclei': Nuclei measurements
        - 'Expand_Nuclei': Expanded nuclei measurements (with Gini, moments, etc.)
        - 'Perinuclear_region': Perinuclear region measurements
        - 'edge_spots': Edge spot measurements (if any)
        - 'Image': Image-level statistics
    """
    # Use masked image for edge_spots if provided, otherwise fall back to miro_image
    edge_spots_miro = masked_miro_image if masked_miro_image is not None else miro_image
    results = {}

    # Image metadata columns
    image_cols = {
        "ImageNumber": metadata.get("ImageNumber", 1),
        "FileName_Hoechst": metadata.get("FileName_Hoechst", ""),
        "FileName_MIRO160mer": metadata.get("FileName_MIRO160mer", ""),
        "PathName_Hoechst": metadata.get("PathName_Hoechst", ""),
        "PathName_MIRO160mer": metadata.get("PathName_MIRO160mer", ""),
    }

    # 1. Nuclei measurements (Hoechst intensity)
    nuclei_df = measure_all_object_properties(
        nuclei_labels,
        hoechst_image,
        channel_name="Hoechst",
        include_edge_intensity=False,
        include_mass_displacement=False,
        include_gini=False,
        include_moments=False,
    )

    # Add metadata
    for col, val in image_cols.items():
        nuclei_df[col] = val

    # Add Number_Object_Number column (for compatibility)
    nuclei_df["Number_Object_Number"] = nuclei_df["ObjectNumber"]

    results["Nuclei"] = nuclei_df

    # 2. Expanded nuclei measurements (MIRO160mer with Gini, moments, mass displacement)
    expand_nuclei_df = measure_all_object_properties(
        expand_nuclei_labels,
        miro_image,
        channel_name="MIRO160mer",
        include_edge_intensity=True,
        include_mass_displacement=True,
        include_gini=True,
        include_moments=True,
    )

    # Add metadata
    for col, val in image_cols.items():
        expand_nuclei_df[col] = val

    results["Expand_Nuclei"] = expand_nuclei_df

    # 3. Perinuclear region measurements (MIRO160mer with all metrics)
    perinuclear_df = measure_all_object_properties(
        perinuclear_region_labels,
        miro_image,
        channel_name="MIRO160mer",
        include_edge_intensity=True,
        include_mass_displacement=True,
        include_gini=True,
        include_moments=True,
    )

    # Add metadata
    for col, val in image_cols.items():
        perinuclear_df[col] = val

    results["Perinuclear_region"] = perinuclear_df

    # 4. Edge spots measurements (Module 12: on Subtract_perinucleus_Miro160mer)
    if edge_spots_labels.max() > 0:
        edge_spots_df = measure_all_object_properties(
            edge_spots_labels,
            edge_spots_miro,  # Use masked MIRO per CellProfiler Module 12
            channel_name="Subtract_perinucleus_Miro160mer",
            include_edge_intensity=True,
            include_mass_displacement=True,
            include_gini=False,
            include_moments=False,
        )

        # Add metadata
        for col, val in image_cols.items():
            edge_spots_df[col] = val

        results["edge_spots"] = edge_spots_df
    else:
        # No edge spots detected
        results["edge_spots"] = pd.DataFrame()

    # 5. Image-level measurements
    image_meas = {}
    image_meas.update(ImageMeasurements.measure_image_intensity(hoechst_image))

    # Rename to include channel name
    image_meas_hoechst = {
        k.replace("Image_Intensity", "Image_Intensity_Hoechst"): v for k, v in image_meas.items()
    }

    # MIRO160mer image measurements
    miro_meas = ImageMeasurements.measure_image_intensity(miro_image)
    image_meas_miro = {
        k.replace("Image_Intensity", "Image_Intensity_MIRO160mer"): v for k, v in miro_meas.items()
    }

    # Combine
    image_df = pd.DataFrame([{**image_cols, **image_meas_hoechst, **image_meas_miro}])
    results["Image"] = image_df

    return results
