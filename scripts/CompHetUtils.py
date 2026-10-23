#' retrieve compound heterozygous pairs from silvar output
import os
import warnings
import pandas as pd
from scripts import FileUtils
    
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

def create_comphet_table(silvar_comphet_vcf_f, out_f):
    assert os.path.isfile(silvar_comphet_vcf_f)
    silvar_comphet_vcf_df = FileUtils.read_vcf(silvar_comphet_vcf_f)
    comphet_df = silvar_comphet_vcf_df.apply(lambda row: retrieve_comphet(row["INFO"]),
                                      axis=1)
    
    comphet_df = pd.concat([silvar_comphet_vcf_df[["CHROM", "POS", "REF", "ALT"]],
               comphet_df], axis=1)
    assert len(comphet_df) == len(silvar_comphet_vcf_df), "Inaccurate number of vavriants"
    
    with open(out_f, "w") as out_fh:
        comphet_df.to_csv(out_fh, sep="\t", index=None, header=True)
    assert os.path.isfile(out_f) & (os.path.getsize(out_f) > 0), "Output file not generated"
    
    return comphet_df

