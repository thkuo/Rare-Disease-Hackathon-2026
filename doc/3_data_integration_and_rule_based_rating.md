# 3. Data integration and rule-based rating

Throughout previous stages, the variation consequence and allele frequency data have become available. Overlaps between
the patient's variants and published MVA-associated variants or genes are also identified. Furthermore, the original VCF file includes quality information. In this stage, the informations will be integrated to score variants. 

Steps of this stage are explained below with Python scripts. 

## Loading data

The input data of this step include: bed-formatted overlapped variants (`matched_var_f`), allele frequecy table (`alt_freq_f`), and consequence table (`conseqeunce_tab_f`). The files are row-wise concatenated data of variant- and gene-based outcomes from previous steps.
```py
import os
import warnings
import random
import pandas as pd

os.getcwd()
# input data
matched_var_f = "VG_overlapped.txt"
alt_freq_f = "af_VG_overlapped.txt"
consequence_tab_f = "consequence_VG_overlapped.txt"
# output data
output_f = "draft_integrated_result_v3.tsv"

# Phase 1: load all tables
# layer 1-2: QC and overlap with published MVA-relevant variants
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
           "end", "name"]
matched_var_df.columns = vcf_col + bed_col

# extract the quality table
qual_df = matched_var_df[["CHROM", "POS", "REF", "ALT", "FILTER"]]
count_df = matched_var_df.groupby(by=["CHROM", "POS", "REF", "ALT"]).size().reset_index(name="count")

# layer 3: allele frequnecy
alt_freq_df = pd.DataFrame()
with open(alt_freq_f, "r") as alt_freq_fh:
    alt_freq_df = pd.read_csv(alt_freq_fh,
                                 header = 0,
                                 delimiter="\t",
                                 index_col = None)

# layer 4: variant consequence
with open(consequence_tab_f, "r") as consequence_tab_fh:
    consequence_df = pd.read_csv(consequence_tab_fh, header=0, sep="\t")
consequence_scores_mapping_df = pd.DataFrame({"impact": ["HIGH", "MODERATE", "MODIFIER"],
                                   "consequence_score": [1.0, 0.5, 0.0]})
consequence_df = consequence_df.merge(consequence_scores_mapping_df, on="impact")
```

## Integrate tables

The tables are joined based "CHROM", "POS","REF", and "ALT" columns. 

```py
id_cols = ["CHROM", "POS","REF", "ALT"]
tables_to_be_integrated = {
    "QC": qual_df,
    "count": count_df,
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
```

## Rate the variants

Each variant is score based on the rules below using integrated information:

```py
# determine the final score
integrated_df = integrated_df.fillna(-1)
# # For all variants overlap with any published causal variant;
# # that is, all in the matched table.
# integrated_df["matched_score"] = 1
# the score is no longer binary but the number of overlapped reference loci
integrated_df["matched_score"] = integrated_df["count"]
# 1, if passing the original QC
integrated_df["qual_score"] = integrated_df["FILTER"].apply(
    lambda x: 1 if x == "PASS" else 0)
# 1, if lower than the cutoff of 0.01
# as the variant is expected to be rare for the rare, inheritable disease
integrated_df["freq_score"] = integrated_df["AF"].apply(
    lambda x: 1 if x <= 0.01 else 0)
integrated_df["integrated_score"] = integrated_df.apply(
    lambda row: row["matched_score"] * row["qual_score"] * row["freq_score"] * row["consequence_score"],
    axis=1)
integrated_df = integrated_df.sort_values(by="integrated_score", ascending=False)

with open(output_f, "w") as output_fh:
    integrated_df.to_csv(output_fh, index=False, sep="\t")
assert (os.path.isfile(output_f)) & (os.path.getsize(output_f) > 0)
```
