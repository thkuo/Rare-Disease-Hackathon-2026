import os
import pandas as pd

def merge_scores_and_comphet(integrated_results_f: str, comphet_tab_f: str):
    """Load and combine tables of integrated scores and comp-hets

    Args:
        integrated_results_f (str): filename of the integrated scores
        comphet_tab_f (str): filename of the table of compound heterozygosity

    Returns:
        DataFrame 
    """
    integrated_df = pd.DataFrame()
    with open(integrated_results_f, "r") as integrated_results_fh:
        integrated_df = pd.read_csv(integrated_results_fh, sep="\t", header=0, index_col=False)

    # include compound heterozygosity data
    comphet_df = pd.DataFrame()
    with open(comphet_tab_f, "r") as comphet_tab_fh:
        comphet_df = pd.read_csv(comphet_tab_fh, sep = "\t",
                                 header=0, index_col=None)
    integrated_df = integrated_df.merge(
        comphet_df, how="left", on=["CHROM", "POS", "REF", "ALT"])

    return integrated_df


def drop_redundant_rows(integrated_df: pd.DataFrame, target_metric: str="prob"):
    """Keep only one row for variants forming compound heterozygosity

    Args:
        integrated_df (pd.DataFrame): the table of integrated score and comp-hets.
        target_metric (str, optional):
            The metric for determing the representative
            row for two variants. Defaults to "prob".

    Returns:
        DataFrame
    """
    # First, sort data based on the target metric
    sorted_integrated_df = integrated_df.sort_values(
        target_metric, ascending=False).reset_index(drop=True)
    sorted_integrated_df["ix"] = sorted_integrated_df.index
    # retrieve variants paired with another in the data
    paired_df = sorted_integrated_df[["CHROM2", "POS2", "REF2", "ALT2", "ix"]].rename(
        columns={"CHROM2": "CHROM", "POS2":"POS", "REF2":"REF", "ALT2":"ALT"})
    paired_df = paired_df[paired_df["CHROM"].notna()]
    # ensure datatypes of the two columns for the searching procedure in the next step
    for c in ["CHROM", "POS"]:
        paired_df[c] = paired_df[c].apply(lambda x: str(int(x)))
        sorted_integrated_df[c] = sorted_integrated_df[c].apply(lambda x: str(int(x)))
    # Then, identify the row where the paired variant is
    sorted_integrated_df = sorted_integrated_df.merge(paired_df, on=["CHROM", "POS", "REF", "ALT"],
                                                      how = "left")
    # Now, ix_x and ix_y are row indices of paired variants.
    is_redundant = sorted_integrated_df.apply(
        lambda row: (row["ix_y"] >= 0) & (row["ix_y"] < row["ix_x"]),
        axis=1)
    nonredundant_integrated_df = sorted_integrated_df[is_redundant.apply(lambda x: not x)]
    # remove temporary columns
    nonredundant_integrated_df = nonredundant_integrated_df.drop(
        columns=["ix_x", "ix_y"]).reset_index(drop=True)
    assert len(nonredundant_integrated_df) <= len(integrated_df)
    assert len(nonredundant_integrated_df.columns) == len(integrated_df.columns)

    return nonredundant_integrated_df

def load_submission_template(submission_template_f: str):
    """Load the template for submission of putatively causal variants

    Args:
        submission_template_f (str): filename of the template.

    Returns:
        DataFrame
    """

    # format the output for submission
    submission_template_df = pd.DataFrame()
    with open(submission_template_f, "r") as submission_template_fh:
        submission_template_df = pd.read_csv(submission_template_fh)
    return submission_template_df

def create_submission_ready_table(overwrite_result=False,
                                  out_f: str = "",
                                  nonredundant_integrated_df=pd.DataFrame(),
                                  submission_template_df=pd.DataFrame()
                                  ):
    """Filter discoveries and reformat the data table

    Args:
        overwrite_result (TYPE, optional):
            enable to write results to an exisitng file. Defaults to False.
        out_f (str, optional):
            output filename. Defaults to "".
        nonredundant_integrated_df (TYPE, optional):
            the integrated scores table. Defaults to pd.DataFrame().
        submission_template_df (TYPE, optional):
            the table for submission. Defaults to pd.DataFrame().

    Returns:
        DataFrame
    """

    # check output file first
    if not overwrite_result:
        assert not os.path.isfile(out_f), (
            f"{out_f} exists\nPlease enable overwrite_write "
            "or change the output filename"
            )

    # extract and rename columns
    reformatted_integrated_df = nonredundant_integrated_df[
        ["CHROM", "POS","REF", "ALT",
         "CHROM2", "POS2","REF2", "ALT2",
         "prob"]]
    reformatted_integrated_df = reformatted_integrated_df.rename(
        columns={"CHROM": "chrom_1", "POS": "pos_1",
                 "REF": "ref_1", "ALT": "alt_1",
                 "CHROM2": "chrom_2", "POS2": "pos_2",
                 "REF2": "ref_2", "ALT2": "alt_2", "prob": "epcr"}
        )
    reformatted_integrated_df["chrom_1"] = reformatted_integrated_df["chrom_1"].apply(
        lambda x: "chr"+str(x))
    reformatted_integrated_df["proband_id"] = "PROBAND01"
    reformatted_integrated_df["finding_type"] = "primary"
    # drop the example rows
    submission_df = submission_template_df.iloc[2:, :].merge(
        reformatted_integrated_df,
        how = "outer", left_index=False, right_index=False)
    submission_df = submission_df.sort_values("epcr", ascending=False)
    # exclude those with a probability of zero
    submission_df = submission_df[submission_df["epcr"] > 0]
    submission_df = submission_df.reset_index(drop=True)

    # check the format
    assert len(submission_df) <= 10, "exceeding the max number"
    for col in submission_template_df.columns:
        assert col in submission_df.columns, f"missing column '{col}'"
    print("Checks passed")

    if out_f != "":
        with open(out_f, "w") as out_fh:
            submission_df.to_csv(out_fh, sep=",", index=False)

    return submission_df
