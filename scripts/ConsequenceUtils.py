import warnings
import pandas as pd
from scripts import FileUtils

def search_attribute(info_field: str, target_attribute: str):
    """Find the value of specific attribute in INFO

    Args:
        info_field (str): the full string of INFO
        target_attribute (str): the attribute to be extracted

    Returns:
        the value of the specific attribute
    """

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
    """Find genetic consequence based on the INFO content

    Args:
        info_field (str): the full string of INFO

    Returns:
        Series of the consequence, impact (i.e., severity), and the gene symbol
    """

    ann=search_attribute(info_field, "ANN")
    representative_ann = ann.split(",")[0]
    representative_ann_fields = representative_ann.split("|")
    consequence_so = representative_ann_fields[1]
    impact = representative_ann_fields[2]
    gene_name = representative_ann_fields[3]
    return pd.Series({"consequence":consequence_so,
                      "impact": impact,
                      "gene": gene_name})

def make_consequence_table(vcf_f: str, out_f: str=None):
    """Create a table of genetic consequences for the given variant data from the VCF file

    Args:
        vcf_f (str): VCF filename
        out_f (str, optional): output filename

    Returns:
        DataFrame of the consequences
    """

    vcf_df = FileUtils.read_vcf(vcf_f)
    consequence_df = vcf_df.apply(
        lambda row: retrieve_variant_consequence(row["INFO"]),
        axis=1)
    consequence_df = pd.concat([vcf_df[["CHROM", "POS", "REF", "ALT"]],
                                consequence_df], axis=1)
    assert len(consequence_df) == len(vcf_df), "Inaccurate number of vavriants"

    if out_f is not None:
        print("Saving result...")
        with open(out_f, "w") as out_fh:
            consequence_df.to_csv(out_fh, sep="\t", index=None, header=True)

    return consequence_df
