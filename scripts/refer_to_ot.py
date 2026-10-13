import os
import pandas as pd

os.chdir("scripts")

# Read Open Targets data of ClnVar
ot_d = "../data/OT_MVA_variants"
ot_df = pd.DataFrame() # aggregate gene-wise data into one table
for f in os.listdir(ot_d):
    print(f)
    full_f = os.path.join(ot_d, f)
    gene_symbol = f.replace(".json", "")
    with open(full_f, "r") as ot_fh:
        target_df = pd.read_json(ot_fh)
        target_df["symbol"] = gene_symbol
        assert "diseaseId" in target_df.columns, "unexpected format"
        
    if len(ot_df) == 0:
        ot_df = target_df
    else:
        ot_df = pd.concat([ot_df, target_df], ignore_index=True)
    
assert ot_df.shape[0] > 0, "empty data"

# retain only confident data
# ot_df["reviewStatus"].value_counts()
suspicious_status = [
    "criteria provided, conflicting classifications",
    "no assertion criteria provided"]
filtered_ot_df = ot_df[
    ot_df["reviewStatus"].apply(lambda x: x not in suspicious_status)]
assert len(filtered_ot_df) > 0
# focus on those with certainty of pathogenicity
def filter_cln_sig(sig_classes: list = []):
    # Given a list of clinical significance labels, 
    # check whether any of them is one of the labels of interests
    sel_cln_significance = [
        'pathogenic',
        'likely pathogenic',
        'likely pathogenic, pathogenic']
    out = any([x in sel_cln_significance for x in sig_classes])
    return out
    

filtered_ot_df = filtered_ot_df[
    filtered_ot_df["clinicalSignificances"].apply(
        lambda x: filter_cln_sig(x))]
assert len(filtered_ot_df) > 0
assert filtered_ot_df.shape[0] <= ot_df.shape[0], "unexpected filtration result"

filtered_ot_df["symbol"].value_counts()

# Some variants are associated with more than one gene due to their impacts
# to a disntant gene. For instance, rs2140912326 is located in BUB1B, which can cause
# a stop change in the gene. It is also a upstream variant to PAK6.
var_rsid_occurence = filtered_ot_df["variantRsId"].value_counts()
var_rsid_occurence[var_rsid_occurence > 1]

# print the rsId as patterns to filter the VCF file
out_f = "causal_var.txt"
with open(out_f, "w") as out_fh:
    var_rsid_occurence.index.to_series().to_csv(
        out_fh, header = False, index = False)
assert os.path.isfile(out_f)