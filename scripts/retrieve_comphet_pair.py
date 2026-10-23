#' retrieve compound heterozygous pairs from silvar output
import os
import warnings
import re
import pandas as pd

os.chdir("scripts/")

def read_vcf(silvar_comphet_vcf_f: str):
    matched_vcf_df = pd.DataFrame()
    with open(silvar_comphet_vcf_f, "r") as silvar_comphet_vcf_fh:
        # The default reading function accepts only one character 
        # as the comment sign, so additional processing is needed here
        vcf_content = [l.strip() for l in silvar_comphet_vcf_fh.readlines()
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

def retrieve_comphet(info_field: str):
    comphet=search_attribute(info_field, "slivar_comphet")
    comphet_fields = comphet.split("/")
    # add an suffix to avoid confusion with the partner
    return pd.Series({"CHROM2": comphet_fields[3],
                      "POS2": comphet_fields[4],
                      "REF2": comphet_fields[5],
                      "ALT2": comphet_fields[6]})

silvar_comphet_vcf_f = "compound_het_annotation_clnvar_matched_filtered.vcf"
assert os.path.isfile(silvar_comphet_vcf_f)
silvar_comphet_vcf_df = read_vcf(silvar_comphet_vcf_f)
comphet_df = silvar_comphet_vcf_df.apply(lambda row: retrieve_comphet(row["INFO"]),
                                  axis=1)

comphet_df = pd.concat([silvar_comphet_vcf_df[["CHROM", "POS", "REF", "ALT"]],
           comphet_df], axis=1)
assert len(comphet_df) == len(silvar_comphet_vcf_df), "Inaccurate number of vavriants"

out_f = "comphet_clnvar_matched.txt"
with open(out_f, "w") as out_fh:
    comphet_df.to_csv(out_fh, sep="\t", index=None, header=True)
assert os.path.isfile(out_f) & (os.path.getsize(out_f) > 0), "Output file not generated"