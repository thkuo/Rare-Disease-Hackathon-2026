# Rare Disease, Real Kid: The MVA Hackathon 2026

## Prerequisite
The files `mva_hackathon.json` and `mva_hackathon_genetics.json` are conda environments used in this study. 

## Repository structure
```console
|-data/ external data from the event organizer or public database
|-scripts/ methods for identifying causal variants
|-dev/ prototypes of methods
```

## Method Overview

The rule-based workflow includes 4 stages. 
1. Knowledge base formation

This stage includes gathering of variant-level and gene-level association with MVA. For variant-level information, MVA-associated variants are identified and downloaded from ClinVar. For gene-level data, MVA-associated genes are identified from MedGen. _hyperlink to detailed doc_

2. Bioinformatic meta-analysis

This stage focuses on bioinformatic analysis to aggregate information for identifying putative causal variants. First, the intersect between the patient's variants (that is, the VCF file provided by the organizer) and MVA-associated variants or genes is computed using bedtools. For each overlapped variant (or hit variant), the number of variant-level or gene-level overlapps are counted. Then, genetic consequences (e.g., missense) are annotaed using snpeff with based on human genome GRCh38p14. Furthermore, alternate allele frequencies of overlapped variants are identified using GnomAD GraphicAPI.

In addition, compound heterozugosity is checked among high-quality, impactful, heterozygous variants using slivar. 

3. Data integration and rule-based rating

Meta-data are parsed from outputs of computational tools mentioned above, followed by re-formatting into analysis-friendly tabular data. With the information, 4 scores are computed for each hit variant: (1) The Quality Score, determined by the FILTER value in the VCF file, represents the confidence of this signal. (2) The Association Score is determined by the number of overlapps with MVA-associated variants or genes, which suggests strength of linkage with sophisticated data. (3) The Frequency Score categorizes variants based on whether their rarity (cutoff: 0.01), which provides an epidemiological view on putative molecular causes to the rare disease. (4) The Consequence Score is determined based on severity of the annotated consequence, which indicates the molecular pathogenecity.

The integrated score is the product of the 4 sub-scores. Then, for each hit variant, bootstrapping is conducted to estimate the likelihood _l_ of obtaing a higher integrated score randomly. The representative probability _P_ is calculated by
    _P=1-l_.

Variants with a probability of 0, whcih are usually those with an integrated score of 0, are filtered. For the rest, their sub-scores and the probiblity are then formatted into a table.

4. Data formatting

Compound heterozygosity (comp-het) information is added to the integrated score table. As the integrated table keeps variants forming comp-het in 2 rows, redundanty rows are then dropped after comp-het data is appended to the table.  

## Author
- Tzu-Hao Kuo (thkuo) [GitHub Profile](https://github.com/thkuo)
