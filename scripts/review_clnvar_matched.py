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

vep_f = "../data/2026-10-14-10-02_VEP.tsv.gz"
vep_df = pd.DataFrame()
with open(vep_f, "rb") as vep_fh:
   vep_df = pd.read_csv(vep_fh, compression="gzip", 
                        header = 0, index_col = None, sep = "\t")
   
# expand the predicted effects
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
        match impact_str:
            case "HIGH":
                return 1.0
            case "MODERATE":
                return 0.5
            case _:
                return 0

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
        return vep_score_df[["CHROM", "POS", "REF", "ALT", 
                            "Consequence", "IMPACT", "consequence_score"]]

vep_obj = vep(vep_df)
vep_obj.preprocess()
vep_obj.make_consequence_scores_table()