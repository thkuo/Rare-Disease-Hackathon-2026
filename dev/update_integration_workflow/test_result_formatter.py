import os
from scripts import SubmissionPrep as sp
os.getcwd()
# input files
overwrite_result=True # to protect the old version or not
integrated_results_f = "dev/update_integration_workflow/draft_integrated_result_v3.tsv"
comphet_tab_f = "dev/update_integration_workflow/comphet_VG_overlapped.txt"
submission_template_f = os.path.join(
    "data", "submission_template",
    "track1_submission_template.csv")
# output file
test_out_f = "dev/update_integration_workflow/thkuo-RuleBased_03.csv"


def test_pipeline(integrated_results_f, comphet_tab_f,
                  submission_template_f, overwrite_result, test_out_f):
    integrated_df = sp.merge_scores_and_comphet(integrated_results_f, comphet_tab_f)
    nonredundant_integrated_df=sp.drop_redundant_rows(integrated_df, target_metric="prob")
    submission_template_df = sp.load_submission_template(submission_template_f)
    test_submission_df = sp.create_submission_ready_table(overwrite_result,
                                                 test_out_f,
                                                 nonredundant_integrated_df,
                                                 submission_template_df)
    return test_submission_df

test_out_df = test_pipeline(integrated_results_f,
                            comphet_tab_f,
                            submission_template_f,
                            overwrite_result,
                            test_out_f)

