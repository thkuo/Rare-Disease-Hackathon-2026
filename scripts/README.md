# Environment and Prerequesite
- Ubuntu
22.04.5 LTS

- Conda
26.5.3

- The other dependencies can be installed by setting up the conda environment:
```sh
conda env create -f ../mva_hackackthon.json
conda env create -f ../mva_hackackthon_genetics.json
conda activate mva_hackthon
```

# Method

The rule-based workflow aims to identify causal variants using public information and the precomputed VCF file. The procedures include:

    i. detect overlaps with published causal variants
    First, generate the bed file `clnvar_mva.bed` that contains MVA-associated variants downloaded from ClinVar. The
    search key word in ClinVar was [mosaic variegated aneuploidy](https://www.ncbi.nlm.nih.gov/clinvar/search/?assembly=GRCh38&disease=mosaic%20variegated%20aneuploidy).
```sh
conda activate mva_hackathon
python refer_to_clnvar.py 
conda deactivate
```
    Then, the file can be used to filter the patient's VCF data and identify variants overlapping or identicle with well-studied causal variants. 
```sh
conda activate mva_hackathon
bedtools intersect -wa -wb \
    -a ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
    -b clnvar_mva.bed \
    > clnvar_mva_matched.bed
conda deactivate
```
    ii. annotate further information for the matched variants
        a. add the variant consequence using snpEFF
```sh
zcat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep '#' > tmp.vcf
cut -f -10 clnvar_mva_matched.bed | sed s/\"//g >> tmp.vcf
snpEff ann -Xmx4g -canon -v GRCh38.p14 -o vcf tmp.vcf > annotation_clnvar_matched.vcf
```
Then, consequence and impact attributes are extracted:
```sh
python retrieve_consequence.py
```

        b. check allele frequency by manual review
        The population-wise allele frequency could vary by the sample size of each genetic study. Additionally, manual
        review or curation could help identify potential inconsistencies or errors in the workflow. Therefore, this step
        visits online databases, such as dbSNP and ClinVar, to aggregate allele frequency information for the variants.

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
    iv. identify putative compound heterozygosity

```sh
# filter variants to focus on high-quality, impactful variants
slivar expr --vcf annotation_clnvar_matched.vcf --pass-only \
    --info 'INFO.impactful && FILTER == "PASS" && s.het' | uniq > annotation_clnvar_matched_filtered.vcf
# identify compound het
slivar compound-hets --allow-non-trios --ped tmp.ped.txt \
    -v annotation_clnvar_matched_filtered.vcf \
    > compound_het_annotation_clnvar_matched_filtered.vcf
# extract the compound heterozygosity information and format them into a table    
python retrieve_comphet_pair.py 
```
The file `tmp.ped.txt` was manually made to describe the singleton (that is, only the patient instead of trios)

    v. format the output data
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
