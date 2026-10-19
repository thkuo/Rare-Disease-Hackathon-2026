
# Method
1. Variant-oriented search using causal variants from Open Targets Platform
```sh
python refer_to_ot.py
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class1.txt
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class2.txt
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class3.txt
```
There is no match.

2. Variant-oriented search using causal variants from ClinVar database

    i. detect overlaps with published causal variants
    _NOTE: The method suggests a limitation in detecting novel causal variants_
```sh
python refer_to_clnvar.py # generating clnvar_mva.bed 
bedtools intersect -wa -wb \
    -a ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
    -b clnvar_mva.bed \
    > clnvar_mva_matched.bed
```
    ii. manually review the matched variants using public data (e.g., dbSNP)
    iii. integrate results
```sh
python review_clnvar_matched.py
```

## Scoring method

- qualification score: reliable signal or not
  1, if QUAL == PASS
  0, otherwise
- association score: overlap with published causal variants
  1, if overlapping with any published causal variants
- consequence score: how impactful is the variant
  1, strong severity (e.g., stop codon)
  .5, moderate severity (e.g., upstream)
  0, otherwise
- frequency score: whether the allele is rare
  1, if Alt allele frequency is <= 0.01

# Result
- clnvar_mva.bed: MVA causal variants listed in [ClnVar](https://www.ncbi.nlm.nih.gov/clinvar/search/?assembly=GRCh38&disease=mosaic%20variegated%20aneuploidy)
- manual_review_clnvar_matched_Data.txt: a table of allele frequencies and genetic consequence per variant
- manual_review_clnvar_matched_Ensembl.txt: impact ratings from Ensembls
- manual_review_clnvar_matched_Readme.txt: further details regarding the Data table
- clnvar_mva_matched.bed: patient variants overlapped with the causal variants in the public data
