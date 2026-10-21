slivar expects high-quality, impactful, heterozygous variants. Therefore, the workflow starts with the annotation and
filtration:
Note: `--force` was used to deal with unmatched chromosome names such as M and MT in the gff and fasta files for
mitochondrial data. 
```sh
bcftools csq -s - --ncsq 40 --force \
    -g ../data/GRCh38p14_genes.gff3.bgz -l -f ../data/GRCh38p14_softmasked.fa ./tmp.vcf -O v -o tmp_annotated.vcf

slivar expr --vcf tmp_annotated.vcf --pass-only \
    --info 'variant.FILTER == "PASS" && variant.GT="0/1"' \
    -o tmp_annotated_silvar_filtered.vcf
```
Then, the result is sent to slivar to detect compound heterozugous pairs. Note that there is no trios information.
```sh
slivar compound-hets --allow-non-trios --ped tmp.ped.txt -v tmp_annotated_silvar_filtered.vcf
```
The tmp.ped.txt describes the singleton information:
```txt
fam_unknown     WGS_EX2312012   .       .       .       2
```

The file `tmp_test.vcf` comtains additional variants as controls to assess the pipeline. The result suggest the
same compound heterozygousity, so no false-possitive was found with the quick assessment.

For the snpEFF annotated version, the match variants are filtered first:
```sh
slivar expr --vcf annotation_clnvar_matched.vcf \
    --pass-only --info 'INFO.impactful && FILTER == "PASS" && s.het' \
    -o annotation_clnvar_matched_filtered.vcf
```
Then, the compound heterozygous pairs are detected:
```sh
slivar compound-hets --allow-non-trios --ped tmp.ped.txt -v annotation_clnvar_matched_filtered.vcf > compound_het_annotation_clnvar_matched_filtered.vcf
```

Again, the result still suggest the same pair. 
