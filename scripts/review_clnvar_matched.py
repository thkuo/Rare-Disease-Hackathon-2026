# Review variants matched by ClnVar data
import os
import pandas as pd

os.chdir("scripts/")

matched_var_f = "clnvar_mva_matched.bed"
assert os.path.isfile(matched_var_f)
matched_var_df = pd.DataFrame()
with open(matched_var_f, "r") as matched_var_fh:
    matched_var_df = pd.read_csv(matched_var_fh,
                                 header = None,
                                 delimiter="\t",
                                 index_col = None)
# make header
vcf_col=["CHROM", "POS","ID","REF",
         "ALT","QUAL","FILTER","INFO",
         "FORMAT","WGS_EX2312012"]
bed_col = ["GRCh38Chromosome", "start",
           "end", "Canonical SPDI"]
matched_var_df.columns = vcf_col + bed_col

review_f = "manual_review_clnvar_matched.xlsx"
assert os.path.isfile(review_f)
review_df = pd.DataFrame()
with open(review_f, "rb") as review_fh:
    review_df = pd.read_excel(review_fh, sheet_name="Sheet1")

candidate_causal_var_df = matched_var_df.merge(
    review_df,
    left_on = ["CHROM", "POS","REF", "ALT"],
    right_on = ["CHROM", "POS", "REF", "ALT"])