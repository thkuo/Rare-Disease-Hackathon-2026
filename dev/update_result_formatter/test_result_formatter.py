import os
import pandas as pd
from scripts import SubmissionPrep as sp

def standard_pipeline(integrated_results_f, comphet_tab_f,
                      submission_template_f, overwrite_result, test_out_f):
    # pipeline starts
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
    
    '''
    A pair of compound heterozygous variants is 1 hit, but the pipeline treat each
    variant individually even if they form compound heterozygosity.
    Therefore, we need to exclude redundant rows.
    '''
    # First, sort data based on the target metric
    target_metric = "prob"
    sorted_integrated_df = integrated_df.sort_values(target_metric, ascending=False).reset_index(drop=True)
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
    nonredundant_integrated_df = nonredundant_integrated_df.drop(columns=["ix_x", "ix_y"] # remove temporary columns
                                    ).reset_index(drop=True)
    assert len(nonredundant_integrated_df) <= len(integrated_df)
    assert len(nonredundant_integrated_df.columns) == len(integrated_df.columns)
    
    # format the output for submission
    submission_template_df = pd.DataFrame()
    with open(submission_template_f, "r") as submission_template_fh:
        submission_template_df = pd.read_csv(submission_template_fh)
    
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
    submission_df = submission_template_df.iloc[2:, :].merge( # drop the example rows
        reformatted_integrated_df,
        how = "outer", left_index=False, right_index=False)
    submission_df = submission_df.sort_values("epcr", ascending=False)
    submission_df = submission_df[submission_df["epcr"] > 0]
    submission_df = submission_df.reset_index(drop=True)
    
    # check the format
    assert len(submission_df) <= 10, "exceeding the max number"
    for col in submission_template_df.columns:
        assert col in submission_df.columns, "missing column '{}'".format(col)
    print("Checks passed")
    
    # save the file
    if not overwrite_result:
        assert not os.path.isfile(standard_out_f)
        
    with open(standard_out_f, "w") as out_fh:
        submission_df.to_csv(out_fh, sep=",", index=False)
        
    return submission_df
   
def test_pipeline(integrated_results_f, comphet_tab_f,
                  submission_template_f, overwrite_result, test_out_f):
    integrated_df = sp.merge_scores_and_comphet(integrated_results_f, comphet_tab_f)
    nonredundant_integrated_df=sp.drop_redundant_rows(integrated_df, target_metric="prob")
    submission_template_df = sp.load_submission_template(submission_template_f)
    test_submission_df = sp.create_submission_ready_table(overwrite_result,
                                                 test_out_f,
                                                 nonredundant_integrated_df,
                                                 submission_template_df)
    return test_submission_df

# input files
overwrite_result=True # to protect the old version or not
integrated_results_f = "scripts/draft_integrated_result_v2.tsv"
comphet_tab_f = "scripts/comphet_clnvar_matched.txt"
submission_template_f = os.path.join(
    "data", "submission_template",
    "track1_submission_template.csv")
# output file
standard_out_f = "dev/update_result_formatter/standard_submission_table.csv"
test_out_f = "dev/update_result_formatter/test_submission_table.csv"
standard_out_df = standard_pipeline(integrated_results_f, comphet_tab_f, submission_template_f, overwrite_result, standard_out_f)
test_out_df = test_pipeline(integrated_results_f, comphet_tab_f, submission_template_f, overwrite_result, test_out_f)
assert standard_out_df.equals(test_out_df), "result not reproduced"
print("DONE")
