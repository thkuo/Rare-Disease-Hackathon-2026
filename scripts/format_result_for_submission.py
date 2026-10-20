import os
import pandas as pd

integrated_results_f = "draft_integrated_result.tsv"
integrated_df = pd.DataFrame()
with open(integrated_results_f, "r") as integrated_results_fh:
    integrated_df = pd.read_csv(integrated_results_fh, sep="\t", header=0, index_col=False)

# format the output for submission
submission_template_f = os.path.join(
    "..", "data", "submission_template",
    "track1_submission_template.csv"
    )
submission_template_df = pd.DataFrame()
with open(submission_template_f, "r") as submission_template_fh:
    submission_template_df = pd.read_csv(submission_template_fh)
    
reformatted_integrated_df = integrated_df[["CHROM", "POS","REF", "ALT", "prob"]]
reformatted_integrated_df = reformatted_integrated_df.rename(
    columns={"CHROM": "chrom_1", "POS": "pos_1",
             "REF": "ref_1", "ALT": "alt_1", "prob": "epcr"}
    )
reformatted_integrated_df["chrom_1"] = reformatted_integrated_df["chrom_1"].apply(
    lambda x: "chr"+str(x))
reformatted_integrated_df["proband_id"] = "PROBAND01"
reformatted_integrated_df["finding_type"] = reformatted_integrated_df["epcr"].apply(
    lambda x: "primary" if x == max(reformatted_integrated_df["epcr"]) else "secondary")
submission_df = submission_template_df.iloc[2:, :].merge( # drop the example rows
    reformatted_integrated_df,
    how = "outer", left_index=False, right_index=False)
submission_df = submission_df.sort_values("epcr", ascending=False)
submission_df = submission_df[submission_df["epcr"] > 0]

# check the format
assert len(submission_df) <= 10, "exceeding the max number"
for col in submission_template_df.columns:
    assert col in submission_df.columns, "missing column '{}'".format(col)
print("Checks passed")

# save the file
hugging_face_username = "thkuo"
model_name = "RuleBased"
version = "01"
out_f = "{}-{}_{}.csv".format(hugging_face_username, model_name, version)
with open(out_f, "w") as out_fh:
    submission_df.to_csv(out_fh, sep=",", index=False)
    
