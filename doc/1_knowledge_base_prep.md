# 1 Knolwedge base preparation

## ClinVar data collection and processing

ClinVar is a database providing detailed information regarding genomic variation and phenotypic changes. A list of MVA-associated variants is downloaded from ClinVar with the search keyword: '[mosaic variegated aneuploidy](https://www.ncbi.nlm.nih.gov/clinvar/search/?assembly=GRCh38&disease=mosaic%20variegated%20aneuploidy)'. The manually downloaded table is then processed to generate a bed-formatted file using:

```sh
conda activate mva_hackathon
python refer_to_clnvar.py 
conda deactivate
```

## MedGen data collection and processing

MedGen is a database harboring medical genetics. MVA-associated genes are identified in the webpages below:
- type 1-3: https://www.ncbi.nlm.nih.gov/medgen/1641418
- type 4: https://www.ncbi.nlm.nih.gov/medgen/1824040
- type 7: https://www.ncbi.nlm.nih.gov/medgen/1824057

Then, their sequential coordinates are extracted from the [gff3 file of GRCh38p14 downloaded from Ensembl](https://www.ensembl.org/genome/GCA_000001405.29) with:

```sh
export full_gff3=../../data/GRCh38p14_genes.gff3.bgz
export out_f=marker_genes.gff3
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=TRIP13;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=BUB1B;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=CEP57;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=CENATAC;' >> $out_f
zcat $full_gff3 | awk '$3 == "gene"' | grep 'Name=MAD1L1;' >> $out_f
```

The file in a gff3 format is then converted into a bed format with:

```sh
awk '{print $1"\t"$4-1"\t"$5"\t"$9}' $out_f > marker_genes.bed
```

This bed file is then used in the next step to identify patient's variants in the loci.

_NOTE: coordinates in a gff3 format are 1-based, but those of the starting site in a bed format is 0-based._
