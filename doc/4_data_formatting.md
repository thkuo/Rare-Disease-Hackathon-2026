# 4. Data formatting

Previous stages have prioritized variants. This stage will further include compound heterozygosity information to
prepare the table for submission.

In this block, `integrated_results_f` comes from the previous step, which includes prioritized variants and their scores. The comp-het data (`comphet_tab_f`) is a concatenated table of variant- and gene-based detection of comp-het using the method in Stage 2. 

```py
import os
import pandas as pd
from scripts import SubmissionPrep as sp

# input files
overwrite_result=True # to protect the old version or not
integrated_results_f = "draft_integrated_result.tsv"
comphet_tab_f = "comphet_VG_overlapped.txt"
submission_template_f = os.path.join(
    "data", "submission_template",
    "track1_submission_template.csv")
# output file
out_f = "submission_table.csv"

def formatting_pipeline(integrated_results_f, comphet_tab_f,
                  submission_template_f, overwrite_result, test_out_f):
    integrated_df = sp.merge_scores_and_comphet(integrated_results_f, comphet_tab_f)
    nonredundant_integrated_df=sp.drop_redundant_rows(integrated_df, target_metric="prob")
    submission_template_df = sp.load_submission_template(submission_template_f)
    test_submission_df = sp.create_submission_ready_table(overwrite_result,
                                                 test_out_f,
                                                 nonredundant_integrated_df,
                                                 submission_template_df)
    return test_submission_df

test_out_df = formatting_pipeline(integrated_results_f, comphet_tab_f, submission_template_f, overwrite_result, test_out_f)
print("DONE")
```
