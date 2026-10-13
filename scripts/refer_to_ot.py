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

# classify variants by clinical significances
sig_cls_df = pd.DataFrame({
    "clinicalSignificances":
        ["pathogenic",
         "likely pathogenic, pathogenic",
         "likely pathogenic",
         "uncertain significance",
         "likely benign",
         "benign, likely benign",
         "benign"],
        "class":[1, 1, 2, 2, 2, 3, 3]})
ranked_significances = filtered_ot_df[
    ["variantRsId", "clinicalSignificances"]].explode(
        "clinicalSignificances").merge(
            sig_cls_df, on="clinicalSignificances",
            how = "left").sort_values(by = "class") 
vars_top_sig_df = ranked_significances.drop_duplicates(
    subset = "variantRsId", keep = "first")

sig_classes = vars_top_sig_df['class'].value_counts().index.to_list()
for sig_class in sig_classes:
    out_f = "causal_var_class{}.txt".format(sig_class)
    try:
        with open(out_f, "w") as out_fh:
            sel_variants = vars_top_sig_df[
                vars_top_sig_df["class"] == sig_class]["variantRsId"]
            sel_variants = sel_variants[sel_variants.apply(lambda x: x is not None)]
            sel_variants.to_csv(
                out_fh, header = False, index = False)
        assert os.path.isfile(out_f), "Error in saving {}".format(out_f)
    except AssertionError as e:
        print(e)