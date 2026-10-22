
Gene selection based on MedGen
- type 1-3: https://www.ncbi.nlm.nih.gov/medgen/1641418
- type 4: https://www.ncbi.nlm.nih.gov/medgen/1824040
- type 7: https://www.ncbi.nlm.nih.gov/medgen/1824057
```sh
export full_gff3=../../data/GRCh38p14_genes.gff3.bgz
export out_f=marker_genes.gff3
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=TRIP13;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=BUB1B;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=CEP57;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=CENATAC;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=MAD1L1;' >> $out_f
```
Then, convert the format to bed for further filtration.
_NOTE: the coordinates in gff3 are 1-based, but that of the starting site in bed format is 0-based._
```sh
awk '{print $1"\t"$4-1"\t"$5"\t"$9}' marker_genes.gff3 > marker_genes.bed
```

The same method like the causal genes are applied again:
```sh
zcat ../../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep '#' > tmp.vcf
cut -f -10 marker_genes_hit.bed | sed s/\"//g >> tmp.vcf
snpEff ann -Xmx4g -canon -v GRCh38.p14 -o vcf tmp.vcf > annotation_marker_genes_hit.vcf
```
, followed by compound heterozyosity detection:

```sh
slivar expr --vcf annotation_marker_genes_hit.vcf --pass-only \
    --info 'INFO.impactful && FILTER == "PASS" && s.het' | uniq > annotation_marker_genes_hit_filtered.vcf

slivar compound-hets --allow-non-trios --ped ../../scripts/tmp.ped.txt \
    -v annotation_marker_genes_hit_filtered.vcf \
    > compound_het_annotation_marker_genes_hit_filtered.vcf
```


