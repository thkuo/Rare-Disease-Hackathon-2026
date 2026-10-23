import re
import pandas as pd

def read_vcf(vcf_f: str):
    vcf_df = pd.DataFrame()
    with open(vcf_f, "r") as vcf_fh:
        # The default reading function accepts only one character 
        # as the comment sign, so additional processing is needed here
        vcf_content = [l.strip() for l in vcf_fh.readlines()
                       if re.match("##", l) is None]
        vcf_content_tab = [row_string.split("\t") for row_string in vcf_content]
        vcf_df = pd.DataFrame(vcf_content_tab[1:],
                                      columns=vcf_content_tab[0])
        vcf_df = vcf_df.rename(columns = {"#CHROM":"CHROM"})
        
    return vcf_df