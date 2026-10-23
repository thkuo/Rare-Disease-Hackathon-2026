#' identify variant consequence and impact information from the annotated VCF
import os
import warnings
import re
import pandas as pd
from scripts import ConsequenceUtils

os.getcwd()
def standard_pipeline(vcf_f):
    def read_vcf(matched_vcf_f: str):
        matched_vcf_df = pd.DataFrame()
        with open(matched_vcf_f, "r") as matched_vcf_fh:
            # The default reading function accepts only one character 
            # as the comment sign, so additional processing is needed here
            vcf_content = [l.strip() for l in matched_vcf_fh.readlines()
                           if re.match("##", l) is None]
            vcf_content_tab = [row_string.split("\t") for row_string in vcf_content]
            matched_vcf_df = pd.DataFrame(vcf_content_tab[1:],
                                          columns=vcf_content_tab[0])
            matched_vcf_df = matched_vcf_df.rename(columns = {"#CHROM":"CHROM"})
            
            return matched_vcf_df
        
    def search_attribute(info_field: str, target_attribute: str):
        info_list = info_field.split(";")
        info_dict = dict()
        for att_val_str in info_list:
            split_out = att_val_str.strip().split("=")
            if len(split_out) == 2:
                info_dict[split_out[0]] = split_out[1]
            else:
                warnings.warn(f"Unrecognized attribute format: {att_val_str}")
        if target_attribute not in info_dict:
            print(info_dict.keys())
            warnings.warn("Attribute not found")
        else:
            return info_dict[target_attribute]
    
    def retrieve_variant_consequence(info_field: str):
        ann=search_attribute(info_field, "ANN")
        representative_ann = ann.split(",")[0]
        representative_ann_fields = representative_ann.split("|")
        consequence_so = representative_ann_fields[1]
        impact = representative_ann_fields[2]
        gene_name = representative_ann_fields[3]
        return pd.Series({"consequence":consequence_so,
                          "impact": impact,
                          "gene": gene_name})
    

    matched_vcf_df = read_vcf(vcf_f)
    consequence_df = matched_vcf_df.apply(lambda row: retrieve_variant_consequence(row["INFO"]),
                         axis=1)
    consequence_df = pd.concat([matched_vcf_df[["CHROM", "POS", "REF", "ALT"]],
               consequence_df], axis=1)
    assert len(consequence_df) == len(matched_vcf_df), "Inaccurate number of vavriants"
    return consequence_df

annotation_vcf_f = "scripts/annotation_clnvar_matched.vcf" # output from snpeff
standard_out_df = standard_pipeline(annotation_vcf_f)
test_out_df = ConsequenceUtils.make_consequence_table(annotation_vcf_f)
assert test_out_df.sort_values("POS").equals(standard_out_df.sort_values("POS")), "result not reproduced"
