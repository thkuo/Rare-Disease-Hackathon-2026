import os
from scripts import ConsequenceUtils

os.getcwd()

annotation_vcf_f = "dev/marker_gene_regions/annotation_marker_genes_hit.vcf" # output from snpeff
out_f = os.path.join("dev",
                      "test_consequence_function_with_marker_gene_hits",
                      "consequence_marker_genes_hit.txt")
out_df = ConsequenceUtils.make_consequence_table(annotation_vcf_f, out_f=out_f)
