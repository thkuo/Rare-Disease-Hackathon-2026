# ClinVar data collection and processing

ClinVar is a database providing detailed information regarding genomic variation and phenotypic changes. The raw data of MVA-associated variants are downloaded from ClinVar. The search key word was '[mosaic variegated aneuploidy](https://www.ncbi.nlm.nih.gov/clinvar/search/?assembly=GRCh38&disease=mosaic%20variegated%20aneuploidy)'. The manually downloaded table is then processed to generate a bed-formatted file using:

```sh
conda activate mva_hackathon
python refer_to_clnvar.py 
conda deactivate
```

