# Rare Disease, Real Kid: The MVA Hackathon 2026
---
## Table of Contents
- Background
- Repository Structure
- Method Overview
- Author
---
## Background

[Rare Disease, Real Kid: The MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026) is a hackathon aiming to raise awareness of mosaic variegated aneuploidy (MVA) syndrome. The [MVA Society](https://mvasociety.org/) states that there are <50 cases in the world. Patients of this rare genetic disease have a high percentage of cells that contain an inaccurate number of chromosomes. Patients usually face developmental abonormalities, but symptoms could differ by sub-types of the disease despite the rarity. To supoprt families influenced by the disease, this hackathon invites people to analyze genetic data and clinical information of a real case to identify genetic causes or associations to the disease. This repository contains personal methodology addressing this challenge. The rule-based workflow described by this repository aims to integrate public data regarding the disease to score the patient's genetic features based on their associations with MVA. 

Besides LICENSE, please note that this study does not aim to provide any clincal advices, and the author does not hold any liability for damages resulting directly or indirectly from codes in this repository.

## Repository Structure
```console
./
|- scripts/ # functions used in the workflow identifying causal variants
|- doc/ # documentation of the workflow, including rationales and scripts
|- dev/ # prototypes of methods
```
## Method Overview

### Prerequisite
The files `mva_hackathon.json` and `mva_hackathon_genetics.json` are conda environments used in this study. The tested conda version is 26.5.3. Codes were developed in Ubuntu 22.04.5 LTS.

### Workflow

The rule-based workflow identifies putatively causal variants for the patient with 4 stages. 

1. Knowledge Base Preparation ([doc](doc/1_knowledge_base_prep.md))

This stage includes gathering data regarding variant-level and gene-level association with MVA from public databases. For variant-level information, MVA-associated variants are identified and downloaded from ClinVar. For gene-level data, MVA-associated genes are identified from MedGen. Sequential position of MVA-associated genes are then extracted from human genome GRCh38p14. 

2. Bioinformatic Analysis ([doc](doc/2_bioinformatics.md))

This stage focuses on bioinformatic analyses to aggregate evidence for identifying putative causal variants. First, the intersect between the patient's variants (that is, the VCF file provided by the organizer) and MVA-associated variants or genes is computed using bedtools. For each overlapped variant (or hit variant), the number of variant-level or gene-level overlapps are counted. Then, genetic consequences (e.g., missense) are annotaed using snpeff based on human genome GRCh38p14. Alternate allele frequencies of overlapped variants are identified using GnomAD GraphicQL API.

This workflow also investigates compound heterozygosity (comp-het), which is a pair of influential variants at the same locus inherited from different sources. The overlapped variants are first filtered based on their qualities and genetic impacts using 'slivar expr'. Then, putative comp-hets are investigated among filtered variants of the singleton (i.e., no trios) using 'slivar compound-hets'.

3. Data Integration and Rule-Based Rating ([doc](doc/3_data_integration_and_rule_based_rating.md))

Meta-data are parsed from outputs of computational tools mentioned above, followed by re-formatting into analysis-friendly tabular data. With the information, 4 scores are computed for each hit variant: (1) The Quality Score, determined by the FILTER value in the VCF file, represents the confidence of this signal. (2) The Association Score is determined by the number of overlapps with MVA-associated variants or genes, which suggests strength of linkage with sophisticated data. (3) The Frequency Score categorizes variants based on whether their rarity (cutoff: 0.01), which provides an epidemiological view on putative molecular causes to the rare disease. (4) The Consequence Score is determined based on severity of the annotated consequence, which indicates the molecular pathogenecity.

The integrated score is the product of the 4 sub-scores. Then, for each hit variant, bootstrapping is conducted to estimate the likelihood _l_ of obtaing a higher integrated score randomly. The representative probability _P_ is calculated as _P=1-l_.

The sub-scores and probiblity are then formatted into a table.


4. Submission Data Formatting ([doc](doc/4_data_formatting.md))

Compound heterozygosity (comp-het) information is added to the table in this stage. For variants in compound heterozygosity, the representative row will be determined based on the probability: the one with a higher probability is kept as the representative, and the other (that is, the one with a lower probability) is consequently excluded from the table to remove redundancy in the data. Then, this stage sets a probability cutoff of 0 to keep only putative candidates as 'primary' discoveries. 

The stage arranges columns based on the submission template and save the table in a csv format for submission.

## Author

Tzu-Hao Kuo, PhD (Taipei, Taiwan)
- Hugging Face: @thkuo
- LinkedIn: [Tzu-Hao Kuo](https://www.linkedin.com/in/th-kuo) 
