# Review variants matched by ClnVar data
import os
import warnings
import pandas as pd
import re

os.chdir("scripts/")

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
review_f = "manual_review_clnvar_matched.xlsx"
assert os.path.isfile(review_f)
review_df = pd.DataFrame()
impact_rating_df = pd.DataFrame()
with open(review_f, "rb") as review_fh:
    review_df = pd.read_excel(review_fh, sheet_name="Sheet1")
    impact_rating_df = pd.read_excel(review_fh,
                                     sheet_name="Ensembl_impact_rating")
alt_freq_df = review_df[["CHROM", "POS","REF", "ALT", 'Alt Frequency']]
alt_freq_df = alt_freq_df.fillna(value = {"Alt Frequency": -1})

# layer 4: variant consequence
vep_f = "../data/2026-10-14-10-02_VEP.tsv.gz"
vep_df = pd.DataFrame()
with open(vep_f, "rb") as vep_fh:
   vep_df = pd.read_csv(vep_fh, compression="gzip", 
                        header = 0, index_col = None, sep = "\t")
   
# expand the predicted effects
def assign_consequence_score(impact_str: str=""):
    match impact_str:
        case "HIGH":
            return 1.0
        case "MODERATE":
            return 0.5
        case _:
            return 0
        
class vep:
    def __init__(self, vep_df):
        self.raw_vep_df = vep_df
        self.preprocessed_vep_df = pd.DataFrame()

    def preprocess(self):
        assert self.raw_vep_df is not None
        vep_df = self.raw_vep_df
        # seperate chromosome and position columns
        vep_reformatted_loc_df = vep_df.apply(
            lambda row: self.reformat_location(row), axis=1)
        vep_df = pd.concat([vep_df, vep_reformatted_loc_df],
                           axis = 1)

        self.preprocessed_vep_df = vep_df
        print("VEP table preprocessed")

    def assign_consequence_score(self, impact_str: str=""):
        return assign_consequence_score(impact_str)            

    def reformat_location(self, row):
        loc_elements = row["Location"].split(":")
        try:
            loc_elements = [int(l) for l in loc_elements]
            return(pd.Series({"CHROM":loc_elements[0], "POS":loc_elements[1]}))
        except ValueError as e:
            print("unrecognized locational label")
            print(e)
            return(pd.Series({"CHROM":-1, "POS":-1}))

    def make_consequence_scores_table(self):
        assert self.preprocessed_vep_df is not None
        vep_df = self.preprocessed_vep_df
        # assign impact scores
        vep_df["consequence_score"] = vep_df.apply(lambda row: self.assign_consequence_score(row["IMPACT"]), axis = 1)
        vep_score_df = vep_df.sort_values(
            by="consequence_score", ascending=False).drop_duplicates(
                subset=["CHROM", "POS", "Ref", "Allele"], keep = "first")
        
        print(vep_score_df.columns)
        vep_score_df = vep_score_df.rename(
            columns={"Ref": "REF", "Allele": "ALT"})
        return vep_score_df[["CHROM", "POS", "REF", "ALT", "SYMBOL",
                            "Consequence", "IMPACT", "consequence_score"]]

vep_obj = vep(vep_df)
vep_obj.preprocess()
vep_consequence_df = vep_obj.make_consequence_scores_table()
vep_consequence_df["Gene Consequence Pair"] = vep_consequence_df.apply(
    lambda row: row["SYMBOL"] + " : " + row["Consequence"], axis = 1)


def identify_representative_impact_of_manually_reviewed_var(
        impact_rating_df,review_df):
    # First, create the mapping table
    impact_mapping_df = impact_rating_df[["Display term", "IMPACT"]]
    impact_mapping_df["consequence"] = impact_mapping_df.apply(
        lambda row: row["Display term"].lower(), axis = 1)
    impact_mapping_df["consequence_score"] = impact_mapping_df.apply(
        lambda row: assign_consequence_score(row["IMPACT"]),
        axis = 1)
    # process those from manual review 
    sub_review_df = review_df[
        ["CHROM", "POS","REF", "ALT", "Gene and Consequence"]][
            review_df["Gene and Consequence"].isna().apply( # exclude the missing ones
                lambda x: not x)]
    sub_review_df["Gene Consequence Pair"] = sub_review_df.apply(
        lambda row: row["Gene and Consequence"].split("|"), axis = 1
        )
    # map the scores
    exploded_sub_review_df = sub_review_df.explode("Gene Consequence Pair")
    exploded_sub_review_df["consequence"] = exploded_sub_review_df.apply(
        lambda row: re.search(": ([\w ]+)", row["Gene Consequence Pair"]).group(1).lower(),
        axis = 1)
    exploded_sub_review_df = exploded_sub_review_df.merge(
        impact_mapping_df,
        left_on = "consequence",
        right_on = "consequence")
    #  determine the representative by finding the most impactful one
    review_consequence_df = exploded_sub_review_df.sort_values(
        by="consequence_score", ascending=False).drop_duplicates(
            subset=["CHROM", "POS", "REF", "ALT"], keep = "first")
            
    return review_consequence_df   
   

review_consequence_df = identify_representative_impact_of_manually_reviewed_var(
    impact_rating_df, review_df)     
target_columns = ['CHROM', 'POS', 'REF', 'ALT',
                  'Gene Consequence Pair','IMPACT',
                  'consequence_score']
consequence_df = pd.concat(
    [vep_consequence_df[target_columns],
    review_consequence_df[target_columns]], ignore_index=True)

# Integrate the tables
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
    warnings.warn("Table {}: {}".format(key, "duplicate found"))
    cleaned_tables_to_be_integrated[key] = non_redundant_df
    try:
        id_col_presence = all([x in tables_to_be_integrated[key].columns for x in id_cols])
        assert id_col_presence, "missing ID columns"
        
    except AssertionError as e:
        id_col_check = False
        warnings.warn("Table {}: {}".format(key, e))
    
if id_col_check:
    print("Table check passed")

integrated_df = pd.DataFrame()
for key in cleaned_tables_to_be_integrated:
    if len(integrated_df) == 0:
        integrated_df = cleaned_tables_to_be_integrated[key]
    else:
        integrated_df = integrated_df.merge(
            cleaned_tables_to_be_integrated[key],
            how = "outer", on=id_cols)
    print(integrated_df.shape)

# determine the final score
integrated_df = integrated_df.fillna(-1)
integrated_df["qual_score"] = integrated_df["FILTER"].apply(
    lambda x: 1 if x == "PASS" else 0)
integrated_df["freq_score"] = integrated_df["Alt Frequency"].apply(
    lambda x: 1 if x <= 0.01 else 0)
integrated_df["integrated_score"] = integrated_df.apply(
    lambda row: row["qual_score"] * row["freq_score"] * row["consequence_score"],
    axis=1)
integrated_df.sort_values(by="integrated_score", ascending=False)