# 2. Bioinformatic meta-analysis

Stage 1 has collected coordinates of published association between the disease and genetic features (i.e., variant or gene). Consequently, this stage deploys bioinformatic tools to map prior knowledge to the patient's data to identify putative causes of MVA. The workflow below includes bash and Python scripts relevant to the analysis. 

## Genetic overlap counting

First, the bed file from Stage 1 enables filtering the patient's VCF file to identify variants overlapping with published MVA-associated variants or genes. The method using [bedtools](https://bedtools.readthedocs.io/en/stable/) is shown below:

```sh
conda activate mva_hackathon
bedtools intersect -wa -wb \
    -a ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
    -b [bed file from Stage 1] \
    > [variant|gene]_overlapped.bed
conda deactivate
```

## Consequence analysis
As molecular causes of inheritable rare diseases are usually impactful but rare, the next steps will focus on genetic
consequences and population observations to add layers to the causal variant identification. 

The workflow described below annotates consequence information for the overlapped variants (that is, those identified from the above step) using [snpEFF](https://pcingola.github.io/SnpEff/) based on [GRCh38p14](https://www.ensembl.org/genome/GCA_000001405.29):

```sh
conda activate mva_hackathon_genetics
zcat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep '#' > tmp.vcf
cut -f -10 [variant|gene]_overlapped.bed| sed s/\"//g >> tmp.vcf
snpEff ann -Xmx4g -canon -v GRCh38.p14 -o vcf tmp.vcf > annotation_[variant|gene]_overlapped.vcf
conda deactivate
```
Please note that the method focuses on canonical sequences. Investigation into isoforms can be the next direction of causal variation identification. 

## Allele frequency

Allele frequencies are retrieved from [GnomAD](https://gnomad.broadinstitute.org/), which aggregates large-scaled sequencing projects, using the GraphicQL API. Below is the Python method to fetch the database with a VCF file:

```py
import os
from scripts import GnomadUtils

vcf_f = "annotation_[variant|gene]_overlapped.vcf"
vcf_f = os.path.join(out_dir, vcf_f)
vcf_df = GnomadUtils.make_af_table(vcf_f=vcf_f, out_f=vcf_f)
vcf_df = vcf_df.sort_values("AF")
assert os.path.isfile(vcf_f) & (os.path.getsize(vcf_f) > 0)
assert vcf_df.equals(standard_out_df)
print("DONE")   
```

For details about the query, please check scripts/GnomadUtils. 

## Compound heterozygosity identification

Compound heterozygosity means paired variants located at the same genetic region of different (opposite) chromosomes. Identifying the type of variants typically involves trios data (that is, parents and the patient) to track the source of specific genetic features. Alternatively, long-read sequencing technique is also conducive, as the data help to determine whether two variants are located at either the same sequence or opposite chromosomes. However, the dataset does not include genetic data of parents, and the sequencing platform was Illumina NovoSeq (please refer to the [thread](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/25) in the community). Therefore, the method below uses [slivar](https://github.com/brentp/slivar) with population genetic data and the singleton information to determine comp-hets instead.

```sh
conda activate mva_hackathon_genetics
export vcf_f=annotation_[variant|gene]_overlapped.vcf 
export filtered_vcf_f=annotation_[variant|gene]_overlapped_filtered.vcf 
export comphet_vcf_f=compound_het_[variant|gene]_overlapped_filtered.vcf 
bcftools view -i 'GT="het"' $vcf_f -O v | \
    slivar expr --vcf -  --pass-only -g ../../data/slivar/gnomad.hg38.v2.zip \
    -g ../../data/slivar/topmed.hg38.dbsnp.151.zip \
    --info 'INFO.impactful && variant.FILTER == "PASS" && INFO.gnomad_popmax_af < 0.01 && INFO.topmed_af < 0.05'| uniq > $filtered_vcf

slivar compound-hets --allow-non-trios --ped tmp.ped.txt \
    -v $filtered_vcf_f > $comphet_vcf_f
conda deactivate
```

The `tmp.ped.txt` describes the singleton as:
```txt
fam_unknown	WGS_EX2312012	.	.	.	2
```

The population genetics data `gnomad.hg38.v2.zip` and `topmed.hg38.dbsnp.151.zip` can be downloaded from [slivar repository](https://github.com/brentp/slivar/releases/tag/v0.3.4).
