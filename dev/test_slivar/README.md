
Adding population data to slivar usage for filtering variants:

```sh
conda activate mva_hackathon_genetics
export vcf_f=../marker_gene_regions/annotation_marker_genes_hit.vcf 
bcftools view -i 'GT="het"' $vcf_f -O v | \
    slivar expr --vcf -  --pass-only -g ../../data/slivar/gnomad.hg38.v2.zip \
    -g ../../data/slivar/topmed.hg38.dbsnp.151.zip \
    --info 'INFO.impactful && variant.FILTER == "PASS" && INFO.gnomad_popmax_af < 0.01 && INFO.topmed_af < 0.05'| uniq > annotation_marker_genes_hit_filtered.vcf
```

```sh
slivar compound-hets --allow-non-trios --ped ../../scripts/tmp.ped.txt \
    -v annotation_marker_genes_hit_filtered.vcf \
    > compound_het_annotation_marker_genes_hit_filtered.vcf
```
