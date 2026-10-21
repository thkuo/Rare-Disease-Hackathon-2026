# Review variants matched by ClnVar data
import os
import warnings
import random
import pandas as pd

os.chdir("scripts/")

# Phase 1: load all tables
# layer 1-2: QC and overlap with published MVA-relevant variants
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

# extract the quality table
qual_df = matched_var_df[["CHROM", "POS", "REF", "ALT", "FILTER"]]

# layer 3: allele frequnecy
# load the manual review table
review_f = "manual_review_clnvar_matched_Data.txt"
assert os.path.isfile(review_f)
review_df = pd.DataFrame()
with open(review_f, "rb") as review_fh:
    review_df = pd.read_table(review_fh, sep = "\t")
# the impact rating from Ensembl
impact_rating_f = "manual_review_clnvar_matched_Ensembl_impact_rating.txt"
assert os.path.isfile(impact_rating_f)
impact_rating_df = pd.DataFrame()
with open(impact_rating_f, "rb") as impact_rating_fh:
    impact_rating_df = pd.read_table(impact_rating_fh)

alt_freq_df = review_df[["CHROM", "POS","REF", "ALT", 'Alt Frequency']]
alt_freq_df = alt_freq_df.fillna(value = {"Alt Frequency": -1})

# layer 4: variant consequence
consequence_tab_f = "consequence_clnvar_matched.txt"
with open(consequence_tab_f, "r") as consequence_tab_fh:
    consequence_df = pd.read_csv(consequence_tab_fh, header=0, sep="\t")
consequence_scores_mapping_df = pd.DataFrame({"impact": ["HIGH", "MODERATE", "MODIFIER"],
                                   "consequence_score": [1.0, 0.5, 0.0]})
consequence_df = consequence_df.merge(consequence_scores_mapping_df, on="impact")

# Pase 2: Integrate the tables
id_cols = ["CHROM", "POS","REF", "ALT"]
tables_to_be_integrated = {
    "QC": qual_df,
    "alt_freq": alt_freq_df,
    "consequence": consequence_df}
# ensure id columns
cleaned_tables_to_be_integrated = {}
id_col_check = True
for key in tables_to_be_integrated:
    non_redundant_df = tables_to_be_integrated[key].drop_duplicates(
        subset = id_cols)
    # ensure non-redundant data
    all_unique = len(tables_to_be_integrated[key]) == len(non_redundant_df)
    if not all_unique:
        warnings.warn("Table {}: {}".format(key, "duplicate found"))
    cleaned_tables_to_be_integrated[key] = non_redundant_df
    try:
        id_col_presence = all(
            [x in tables_to_be_integrated[key].columns
             for x in id_cols])
        assert id_col_presence, "missing ID columns"

    except AssertionError as e:
        id_col_check = False
        warnings.warn("Table {}: {}".format(key, e))

if id_col_check:
    print("Tables evaluation passed")

integrated_df = pd.DataFrame()
for key in cleaned_tables_to_be_integrated:
    if len(integrated_df) == 0:
        integrated_df = cleaned_tables_to_be_integrated[key]
    else:
        integrated_df = integrated_df.merge(
            cleaned_tables_to_be_integrated[key],
            how = "outer", on=id_cols)
    print(integrated_df.shape)

# Phase 3: Rate the variants
# determine the final score
integrated_df = integrated_df.fillna(-1)
# For all variants overlap with any published causal variant;
# that is, all in the matched table.
integrated_df["matched_score"] = 1
# 1, if passing the original QC
integrated_df["qual_score"] = integrated_df["FILTER"].apply(
    lambda x: 1 if x == "PASS" else 0)
# 1, if lower than the cutoff of 0.01
# as the variant is expected to be rare for the rare, inheritable disease
integrated_df["freq_score"] = integrated_df["Alt Frequency"].apply(
    lambda x: 1 if x <= 0.01 else 0)
integrated_df["integrated_score"] = integrated_df.apply(
    lambda row: row["qual_score"] * row["freq_score"] * row["consequence_score"],
    axis=1)
integrated_df = integrated_df.sort_values(by="integrated_score", ascending=False)

# count p-value by randomization
random.seed(1234)
n = 10000
sample = random.choices(
    integrated_df["integrated_score"],
    k=n)
# The probability obtaining a higher score randomly
integrated_df["p-val"] = integrated_df.apply(
    lambda row: sum([float(s) >= row["integrated_score"] for s in sample])/n,
    axis = 1)
integrated_df["prob"] = integrated_df.apply(
    lambda row: 1 - row["p-val"], axis=1)

output_f = "draft_integrated_result_v2.tsv"
with open(output_f, "w") as output_fh:
    integrated_df.to_csv(output_fh, index=False, sep="\t")
assert (os.path.isfile(output_f)) & (os.path.getsize(output_f) > 0)