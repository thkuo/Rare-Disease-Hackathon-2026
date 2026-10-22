
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
