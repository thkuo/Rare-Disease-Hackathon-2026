from scripts import CompHetUtils as ch

silvar_comphet_vcf_f = "dev/marker_gene_regions/compound_het_annotation_marker_genes_hit_filtered.vcf"
test_out_f = "dev/refactor_comphet_data_parser/test_comphet_marker_genes_hit.txt"

test_df = ch.create_comphet_table(silvar_comphet_vcf_f, test_out_f)