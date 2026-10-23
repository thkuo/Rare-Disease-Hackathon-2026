# 2. Bioinformatic meta-analysis

Stge 1 has collected coordinates of published association between the disease and genetic features (i.e., variant or gene). In this stage, bioinformatic tools are deployed to map prior knowledge to the patient's data to identify putative causes to MVA. The workflow below includes bash and Python scripts relevant to the analysis. 


First, the bed file from Stage 1 enables filtering the patient's VCF file to identify variants overlapping with published causal variants. The method is shown below:

```sh
conda activate mva_hackathon
bedtools intersect -wa -wb \
    -a ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
    -b [bed file from Stage 1] \
    > [variant|gene]_overlapped.bed
conda deactivate
```

As molecular causes to inheritable rare diseases are usually impactful but rare, the next steps will focus on genetic
consequences and population observations to add layers to the causal variant identification. 
The workflow annotates consequence information for the overlapped variants (that is, those identified from the above
step) using snpEFF based on GRCh38p14:

```sh
zcat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep '#' > tmp.vcf
cut -f -10 [variant|gene]_overlapped.bed| sed s/\"//g >> tmp.vcf
snpEff ann -Xmx4g -canon -v GRCh38.p14 -o vcf tmp.vcf > annotation_[variant|gene]_overlapped.vcf
```
Please note that the method focuses on canonical sequences. Investigation into isoforms can be the next direction of causal variation identification. 

Allele frequencies are retrieved from GnomAD, which aggregates large-scaled sequencing projects, using the official GraphicAPI. Below is the Python method to fetch the database with a VCF file:

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
