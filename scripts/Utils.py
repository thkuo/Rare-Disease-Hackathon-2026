import re
import pandas as pd
 
def assign_consequence_score(impact_str: str=""):
     match impact_str:
         case "HIGH":
             return 1.0
         case "MODERATE":
             return 0.5
         case _:
             return 0
  
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