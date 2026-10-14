import os
import re
import pandas as pd

os.chdir("scripts")

# Read data of ClnVar
var_f = "../data/clinvar_download_20261013_040408.tsv"
var_df = pd.DataFrame() # aggregate gene-wise data into one table
with open(var_f, "r") as var_fh:
    var_df = pd.read_csv(var_fh, sep="\t", header=0, index_col = False)
assert var_df.shape[0] > 0

# As the provided VCF was based on GRCh38, those absent from the version are useless
var_df = var_df[(var_df["GRCh38Location"].isna()).apply(lambda x: not x)]

# 
var_df["Germline classification"].value_counts()
var_df["Germline review status"].value_counts()
status_occ = var_df.value_counts(["Germline classification", "Germline review status"])

def reformat_coordinates(location_str: str=""):
    print(location_str)
    try:
        start_loc = int(location_str) - 1 # 0-based start site in BED format
        end_loc = int(location_str)
    except ValueError:
        search_out = re.search("([0-9]+) \- ([0-9]+)", location_str)
        start_loc = int(search_out.group(1)) - 1
        end_loc = int(search_out.group(2))
    
    return pd.Series({"start": start_loc, "end": end_loc})

        
sub_var_df = var_df[["GRCh38Chromosome", "GRCh38Location",
                     "Canonical SPDI",
                     "Germline classification", "Germline review status"]]
sub_var_df[["start", "end"]] = sub_var_df.apply(
    lambda x: reformat_coordinates(x["GRCh38Location"]), axis = 1)
sub_var_df["GRCh38Chromosome"] = sub_var_df["GRCh38Chromosome"].apply(lambda x: int(x))
sub_var_df = sub_var_df[["GRCh38Chromosome", "start", "end", "Canonical SPDI"]] # BED format explained: https://bedtools.readthedocs.io/en/latest/content/general-usage.html

out_bed_f = "clnvar_mva.bed"
with open(out_bed_f, "w") as out_bed_fh:
    sub_var_df.to_csv(out_bed_fh,sep = "\t",
                      header=False, index=False)