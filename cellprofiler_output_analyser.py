"""
Given a folder with the following structure:

INPUT_FOLDER
    - subfolder1
      - All_measurements.csv
      - Expand_nuclei.csv
      - Perinuclear_region.csv

In which either wellnumber and XY vary, or wellnumber, XY, and T vary, do the following:

- from All_measurements, extract the number of edge spots as a fraction of the number of nuclei per XY.
  Generate intermediate with Well, XY, T, and fraction, one row per (well, XY, T) combination.
  Then:

  If T doesn't vary, generate a file with:

    |  Well_1  | Well_2 | Well_3   ....
    | -------------------------
XY1 |  fraction | fraction | fraction
XY2 ....

 If T varies, generate a file for each well number with:

     |  T1      |     T2   |    T3 ....
     | -------------------------
XY1  | fraction | fraction | fraction
XY2  ....

Then average over XYs and plot:

    |    Well_1      |     Well_2    | Well_3   ....
    | -------------------------
T1  |  mean_fraction | mean_fraction | mean_fraction
...

- from Expand_nuclei, extract the mass displacement of the 60mer/mito for each cell. As this is a per-cell measure,
  we don't really care about XY, but track anyway. Generate intermediate with Well, XY, T, and
  mass displacement, where there are N rows per (well, XY, T) combination, where N is the number
  of cells in that field of view. Then:

  If T doesn't vary, generate a file with:

    |   Well_1  |  Well_2   | Well_3   ....
    | XY1 | XY2 | XY1 | XY2 | XY1 | XY2
    | -------------------------
    |  mass | mass | mass | mass | mass | mass
    |  mass | mass | mass | mass | mass | mass
    ...


And
   |  Well_1 |  Well_2   | Well_3   ....
   | -------------------------
   | mass    |   mass    |   mass

with Well_1 contains all XY1 values, stacked.


NB the above is not square since each FoV may have a different number of cells. Fill with NaNs.

If T varies, generate a file for each well number with:

    | T1        |    T2     |    T3 ....
    | XY1 | XY2 | XY1 | XY2 | XY1 | XY2
    | ------------
    |  mass | mass | mass | mass | mass | mass
    ....

Plus the stacked variant as above.

Then average over values and generate:

    | Well_1 | Well_2 | Well_3 ...

 T1 | mean_mass | mean_mass | mean_mass
 ...

- from Perinuclear_region, extract the CoV for each cell, and then generate files as for the mass
  displacement.


Script parameters:
    - input_folder: path to the folder containing the input files
    - output_folder: path to the folder where the output files will be written
    - t_varies: boolean, whether the T parameter varies or not
    - plot: boolean, whether to plot the output.


"""
import glob
import logging
import os
import re

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# ----------- CHANGE HERE ---------------
INPUT_FOLDER = "input_folder"  # the path to all the input files

INPUT_SUBFOLDERS = ["batch17"]
OUTPUT_FOLDER = "output_folder"  # output will be saved to OUTPUT_FOLDER/INPUT_SUBFOLDER
EDGE_SPOT_FILE = "All_measurements.csv"
MASS_DISPLACEMENT_FILE = "Expand_Nuclei.csv"
COV_FILE = "Perinuclear_region.csv"
# Sometimes the mito or 60mer is not present. If not, comment out the line:
MASS_DISPLACEMENT_COLS = {
    # "mass_displacement_mito": "Intensity_MassDisplacement_mito",
    "mass_displacement_60mer": [
        "Intensity_MassDisplacement_MIRO160mer",
        "Intensity_MassDisplacement_MIRO160mer_rescaled",
    ],
    # "mass_displacement_60mer": "Intensity_MassDisplacement_pex",
}
LOGGING_LEVEL = logging.INFO  # logging.INFO or logging.ERROR  normally
PLOT = False
T_VARIES = False

bonus_cols = ["Intensity_MassDisplacement_MIRO160mer", "GINI_Gini_MIRO160mer"]

# ---------------------------------------


def generate_edge_spot_files(
    input_path: str, output_folder: str, t_varies: bool, do_plot=True
):
    """
    Generate the edge_spot_count / nuclei_count for each field of view.

    Assumes two header columns.

    Args:
        - input_path: the path to the input file. Must have the following columns:
        - output_folder: where to output the processed files.
        - t_varies. If we have multiple timepoints and must therefore process the file differently.
            We check this with asserts.
        - do_plot: whether to plot the output.
    """
    logger.info(
        f"Generating edge spot data for {input_path}, output to {output_folder}, t_varies={t_varies}"
    )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_input_df = pd.read_csv(input_path, header=[0, 1])

    processed_df = extract_edgespot_cols(raw_input_df, t_varies)
    intermediate_filepath = os.path.join(output_folder, "edge_spot_fraction_raw.csv")
    logger.info(f"Writing edge spot intermediate to {intermediate_filepath}")
    processed_df.to_csv(intermediate_filepath, index=False)

    if t_varies:
        # File for each well number, T as columns, XY as rows
        for well_number in processed_df.WellNumber.unique():
            well_number_subdf = processed_df[processed_df.WellNumber == well_number]
            output_filename = os.path.join(
                output_folder, f"edge_spot_fraction_over_time_{well_number}.csv"
            )
            save_pivot_table(
                data=well_number_subdf,
                values="edge_spot_fraction",
                index="XY",
                columns="T",
                output_filename=output_filename,
            )

        # average over XYs, normalise and plot:
        average_over_time = (
            processed_df.groupby(["WellNumber", "T"])
            .edge_spot_fraction.mean()
            .reset_index()
        )
        pivot = pd.pivot_table(
            data=average_over_time,
            values="edge_spot_fraction",
            index="T",
            columns="WellNumber",
        )
        pivot = pivot.divide(pivot.iloc[0])
        output_path = os.path.join(
            output_folder, "edge_spot_fraction_mean_over_time_normalised.csv"
        )
        logger.info(
            f"Normalised pivot table has shape {pivot.shape}, writing to {output_path}"
        )
        pivot.to_csv(output_path)
        if do_plot:
            pivot.plot(title="Edge spot fraction over time, normalised to T0")
            plt.show()

    else:
        # single file - well number vs xy.
        output_filename = os.path.join(output_folder, "edge_spot_fraction_static.csv")
        save_pivot_table(
            data=processed_df,
            values="edge_spot_fraction",
            index="XY",
            columns="WellNumber",
            output_filename=output_filename,
        )


def extract_edgespot_cols(
    cellprofiler_df: pd.DataFrame, t_varies: bool
) -> pd.DataFrame:
    """
    Extract the (well_number, xy, t) columns and aggregate over W, XY, T
    to generate the edge spot data.
    Makes a number of assumptions about the struture of the input df.
    """
    logger.info(f"Extracting columns, input df has shape {cellprofiler_df.shape}")
    filename_column = "FileName_Hoechst"  # from which we extract the well number and xy
    nuclei_count_column = "Number_Object_Number"
    edge_spot_count_column = "Number_Object_Number"

    output_df = cellprofiler_df.copy()
    if t_varies:
        output_df["T"] = output_df["Image"][filename_column].apply(extract_timestamp)
    output_df["XY"] = output_df["Image"][filename_column].apply(extract_xy)
    output_df["WellNumber"] = output_df["Image"][filename_column].apply(
        extract_wellnumber
    )
    output_df["nuclei_count"] = output_df["Nuclei"][nuclei_count_column]
    output_df["edge_spot_count"] = output_df["edge_spots"][edge_spot_count_column]
    cols = ["WellNumber", "XY", "nuclei_count", "edge_spot_count"]
    if t_varies:
        cols.append("T")
    output_df = output_df[cols].reset_index(drop=True)
    output_df.columns = output_df.columns.droplevel(
        1
    )  # drop the second level of the multiindex

    aggregate_df = (
        output_df.groupby(
            ["WellNumber", "XY", "T"] if t_varies else ["WellNumber", "XY"]
        )
        .agg({"nuclei_count": "count", "edge_spot_count": "count"})
        .reset_index()
    )
    aggregate_df["edge_spot_fraction"] = (
        aggregate_df["edge_spot_count"] / aggregate_df["nuclei_count"]
    )

    logger.info(
        f"Extracted columns and aggregated, output df has shape {aggregate_df.shape}"
    )
    return aggregate_df


def extract_massdisplacement_cols(cellprofiler_df, t_varies: bool) -> pd.DataFrame:
    """Extract well number, xy, t, mass, displacement columns from the cellprofiler df.
    As above, makes strong assumptions about input shape and labels.

    Christina is a criminal and thus the miro and mito displacments may or may not be there.
    """
    logger.info(f"Extracting columns, input df has shape {cellprofiler_df.shape}")
    # filename_column = "FileName_mito"
    # filename_column = "FileName_pex"
    filename_column = "FileName_MIRO160mer"
    output_df = cellprofiler_df.copy()
    if t_varies:
        output_df["T"] = output_df[filename_column].apply(extract_timestamp)
    output_df["XY"] = output_df[filename_column].apply(extract_xy)
    output_df["WellNumber"] = output_df[filename_column].apply(extract_wellnumber)
    cols = []
    if "mass_displacement_mito" in MASS_DISPLACEMENT_COLS:
        col = MASS_DISPLACEMENT_COLS["mass_displacement_mito"]
        output_df["mass_displacement_mito"] = output_df[col]
        cols.append("mass_displacement_mito")
    if "mass_displacement_60mer" in MASS_DISPLACEMENT_COLS:
        for col in MASS_DISPLACEMENT_COLS["mass_displacement_60mer"]:
            try:
                output_df["mass_displacement_60mer"] = output_df[col]
                break
            except KeyError:
                continue
        else:
            raise ValueError('Could not find a column for "mass_displacement_60mer"')
        cols.append("mass_displacement_60mer")
    cols = ["WellNumber", "XY"] + cols
    if t_varies:
        cols.append("T")
    output_df = output_df[cols].reset_index(drop=True)
    return output_df


def extract_cov_cols(
    cellprofiler_df, t_varies: bool, bonus_cols: list[str] = []
) -> pd.DataFrame:
    """Extract well number, xy, t, and CoV."""
    logger.info(f"Extracting columns, input df has shape {cellprofiler_df.shape}")
    filename_column = "FileName_MIRO160mer"
    std_columns = [
        "Intensity_StdIntensity_MIRO160mer",
        "Intensity_StdIntensity_MIRO160mer_rescaled",
    ]
    mean_columns = [
        "Intensity_MeanIntensity_MIRO160mer",
        "Intensity_MeanIntensity_MIRO160mer_rescaled",
    ]
    output_df = cellprofiler_df.copy()
    if t_varies:
        output_df["T"] = output_df[filename_column].apply(extract_timestamp)
    output_df["XY"] = output_df[filename_column].apply(extract_xy)
    output_df["WellNumber"] = output_df[filename_column].apply(extract_wellnumber)
    for std_col, mean_col in zip(std_columns, mean_columns):
        try:
            output_df["CoV"] = output_df[std_col] / output_df[mean_col]
            break
        except KeyError:
            continue
    else:
        raise ValueError('Could not find a column for "CoV"')
    cols = ["WellNumber", "XY", "CoV"]
    if t_varies:
        cols.append("T")
    cols += bonus_cols
    output_df = output_df[cols].reset_index(drop=True)
    return output_df


def extract_xy(input_string: str) -> int:
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


def extract_wellnumber(input_string: str) -> str:
    return re.search(r"(?<=Well)[A-Z]\d+", input_string).group(0)


def save_pivot_table(data, values, index, columns, output_filename: str):
    pivot = pd.pivot_table(data=data, values=values, index=index, columns=columns)
    logger.info(f"writing pivot table: {output_filename}")
    pivot.to_csv(output_filename)


def generate_mass_displacement_files(
    mass_displacement_file_path, output_folder, t_varies, do_plot=True
):
    """
    Aggregate the mass displacement over all cells.

    Args:
        - input_path: the path to the input file. Must have the following columns:
        - output_folder: where to output the processed files.
        - t_varies. If we have multiple timepoints and must therefore process the file differently.
            We check this with asserts.
        - do_plot: whether to plot the output.
    """
    logger.info(
        f"Generating mass displacement data for {mass_displacement_file_path}, output to {output_folder}, t_varies={t_varies}"
    )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_input_df = pd.read_csv(mass_displacement_file_path)
    processed_df = extract_massdisplacement_cols(raw_input_df, t_varies)
    for displacement_type in MASS_DISPLACEMENT_COLS:
        generate_ragged_df(
            processed_df,
            data_column=displacement_type,
            output_folder=output_folder,
            t_varies=t_varies,
            do_plot=do_plot,
        )


def generate_cov_files(cov_file_path, output_folder, t_varies, do_plot=True):
    """Aggregate the CoV over all cells ( std intensity / mean intensity)

    Highly similar to the mass displacement function.
    """
    logger.info(
        f"Generate CoV data for {cov_file_path}, output to {output_folder}, t_varies={t_varies}"
    )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_input_df = pd.read_csv(cov_file_path)
    processed_df = extract_cov_cols(raw_input_df, t_varies)
    generate_ragged_df(
        processed_df,
        data_column="CoV",
        output_folder=output_folder,
        t_varies=t_varies,
        do_plot=do_plot,
    )


def generate_extra_dispersion_measures(
    cov_file_path, output_folder, t_varies, bonus_cols, do_plot=True
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_input_df = pd.read_csv(cov_file_path)
    processed_df = extract_cov_cols(raw_input_df, t_varies, bonus_cols=bonus_cols)

  
    derived_cols = [  #  insert extra cols here 
    ]

    for bonus_col in bonus_cols + derived_cols:
        generate_ragged_df(
            processed_df,
            data_column=bonus_col,
            output_folder=output_folder,
            t_varies=t_varies,
            do_plot=do_plot,
        )


def generate_ragged_df(
    input_df, data_column, output_folder, t_varies: bool, do_plot=True
):
    """In both mass displacement and CoV we extract some vals and stack them into a ragged df."""
    if t_varies:
        # A file per wellnumber, with T and XY as columns and subcolumns
        subdf = input_df[["WellNumber", "XY", "T", data_column]]
        for well_number in subdf.WellNumber.unique():
            well_number_subdf = subdf[subdf.WellNumber == well_number]
            # highly mysterious pandas esoterica to generate columns for T and XY
            output_filename = os.path.join(
                output_folder, f"{data_column}_time_and_xy_{well_number}.csv"
            )
            output_df = (
                well_number_subdf.drop(columns=["WellNumber"])
                .set_index(["T", "XY"])
                .squeeze()
            )
            output_df = output_df.reset_index(name="value")
            output_df["index"] = output_df.groupby(["T", "XY"]).cumcount()
            output_df = output_df.set_index(["T", "XY", "index"])["value"].unstack(
                level=[0, 1]
            )
            logger.info(
                f"writing {data_column} table with shape {output_df.shape}: {output_filename}"
            )
            output_df.to_csv(output_filename)

            # pivot to generate a column for each T
            stacked_df = well_number_subdf.drop(columns=["WellNumber", "XY"])
            stacked_df["count"] = stacked_df.groupby("T").cumcount()
            stacked_df = stacked_df.pivot(
                index="count", columns="T", values=data_column
            )
            output_filename = os.path.join(
                output_folder, f"{data_column}_stacked_{well_number}.csv"
            )
            logger.info(
                f"writing stacked {data_column} table with shape {stacked_df.shape}: {output_filename}"
            )
            stacked_df.to_csv(output_filename)

        # Average over Well, T, save and plot
        average_over_time = (
            subdf.groupby(["WellNumber", "T"])[data_column].mean().reset_index()
        )
        pivot = pd.pivot_table(
            data=average_over_time,
            values=data_column,
            index="T",
            columns="WellNumber",
        )
        pivot = pivot.divide(pivot.iloc[0])
        output_path = os.path.join(
            output_folder, f"{data_column}_mean_over_time_normalised.csv"
        )
        logger.info(
            f"Normalised pivot table has shape {pivot.shape}, writing to {output_path}"
        )
        pivot.to_csv(output_path)
        if do_plot:
            pivot.plot(title=f"Mean {data_column} over time, normalised to T0")
            plt.show()

    else:
        # Generate a table with WellNumber and XY as columns, and the data column as the vals
        output_filename = os.path.join(output_folder, f"{data_column}_static.csv")
        subdf = input_df[["WellNumber", "XY", data_column]].copy()
        subdf["index"] = subdf.groupby(["WellNumber", "XY"]).cumcount()
        subdf = subdf.pivot(
            index="index", columns=["WellNumber", "XY"], values=data_column
        )
        logger.info(
            f"writing static {data_column} table with shape {subdf.shape}: {output_filename}"
        )
        subdf.to_csv(output_filename)

        # Generate a table with Wellnumber as columns, XY as rows, and the median of the data column as the vals.
        median_df = (
            input_df.groupby(["WellNumber", "XY"])[[data_column]].median().unstack().T
        )
        median_filename = os.path.join(
            output_folder, f"{data_column}_fov_median_static.csv"
        )
        median_df.to_csv(median_filename)

        # Pivot to generate a column for each WellNumber
        stacked_df = input_df[["WellNumber", data_column]].copy()
        stacked_df["count"] = stacked_df.groupby("WellNumber").cumcount()
        stacked_df = stacked_df.pivot(
            index="count", columns="WellNumber", values=data_column
        )
        output_filename = os.path.join(
            output_folder, f"{data_column}_stacked_static.csv"
        )
        logger.info(
            f"writing stacked {data_column} table with shape {stacked_df.shape}: {output_filename}"
        )
        stacked_df.to_csv(output_filename)



if __name__ == "__main__":
    for INPUT_SUBFOLDER in INPUT_SUBFOLDERS:
        input_folders = glob.glob(os.path.join(INPUT_FOLDER, INPUT_SUBFOLDER, "*/"))
        for input_folder in input_folders:
            logger.info(f"Processing {input_folder}")
            output_subfolder = input_folder.replace(INPUT_FOLDER, OUTPUT_FOLDER)
            edge_spot_file_path = os.path.join(input_folder, EDGE_SPOT_FILE)
            generate_edge_spot_files(
                edge_spot_file_path, output_subfolder, T_VARIES, PLOT
            )
            logger.info("\n")
            mass_displacement_file_path = os.path.join(
                input_folder, MASS_DISPLACEMENT_FILE
            )
            generate_mass_displacement_files(
                mass_displacement_file_path, output_subfolder, T_VARIES, PLOT
            )
            logger.info("\n")
            cov_file_path = os.path.join(input_folder, COV_FILE)
            generate_cov_files(cov_file_path, output_subfolder, T_VARIES, PLOT)
            generate_extra_dispersion_measures(
                cov_file_path,
                output_subfolder,
                T_VARIES,
                bonus_cols=bonus_cols,
                do_plot=PLOT,
            )
