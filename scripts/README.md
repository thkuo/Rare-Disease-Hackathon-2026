
# Method
1. Variant-oriented search using causal variants from Open Targets Platform
```sh
python refer_to_ot.py
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class1.txt
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class2.txt
cat ../data/WGS_EX2312012_HGWCNDSX7.vcf.gz | \
grep -w -i -f ./causal_var_class3.txt
```
There is no match.


