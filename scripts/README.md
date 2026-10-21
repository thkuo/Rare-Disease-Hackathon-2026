# Environment and Prerequesite
- Ubuntu
22.04.5 LTS

- Conda
26.5.3

- The other dependencies can be installed by setting up the conda environment:
```sh
conda env create -f ../mva_hackackthon.json
conda activate mva_hackthon
```

# Method

The rule-based workflow aims to identify causal variants using public information and the precomputed VCF file. The procedures include:

    i. detect overlaps with published causal variants
```sh
python refer_to_clnvar.py # generating clnvar_mva.bed 
bedtools intersect -wa -wb \
    -a ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
    -b clnvar_mva.bed \
    > clnvar_mva_matched.bed
```
    ii. annotate the variant consequence using snpEFF
```sh
zcat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep '#' > tmp.vcf
cut -f -10 clnvar_mva_matched.bed >> tmp.vcf
snpEff ann -Xmx4g -canon -v GRCh38.p14 -o vcf tmp.vcf > annotation_clnvar_matched.vcf
```

    iii. integrate results
```sh
python review_clnvar_matched.py
```
    The step involves in scoring of the candidate variants. The rule-based scoring procedure include:
```console
- qualification score: reliable signal or not
  1, if QUAL == PASS
  0, otherwise
- frequency score: whether the allele is rare
  1, if Alt allele frequency is <= 0.01
- association score: overlap with published causal variants
  1, if overlapping with any published causal variants
- consequence score: how impactful is the variant
  1, strong severity (e.g., stop codon)
  .5, moderate severity (e.g., upstream)
  0, otherwise
```
    iv. format the output data
```sh
python format_result_for_submission.py
```

# Result
- clnvar_mva.bed: MVA causal variants listed in [ClnVar](https://www.ncbi.nlm.nih.gov/clinvar/search/?assembly=GRCh38&disease=mosaic%20variegated%20aneuploidy)
- manual_review_clnvar_matched_Data.txt: a table of allele frequencies and genetic consequence per variant
- manual_review_clnvar_matched_Ensembl.txt: impact ratings from Ensembls
- manual_review_clnvar_matched_Readme.txt: further details regarding the Data table
- clnvar_mva_matched.bed: patient variants overlapped with the causal variants in the public data
- draft_integrated_result.tsv: the integrated results with sub-scores
- thkuo-RuleBased_01.csv: final table for the submission
