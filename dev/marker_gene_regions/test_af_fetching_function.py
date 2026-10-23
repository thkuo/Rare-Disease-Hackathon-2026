import os
import requests
import warnings
import time
import re
import pandas as pd
from tqdm import tqdm
from tempfile import TemporaryDirectory, mkdtemp

os.getcwd()
from scripts import GnomadUtils

tqdm.pandas()
def standard_pipeline(vcf_f): 
    def read_vcf(marker_hit_vcf_f: str):
        marker_hit_vcf_df = pd.DataFrame()
        with open(marker_hit_vcf_f, "r") as marker_hit_vcf_fh:
            # The default reading function accepts only one character 
            # as the comment sign, so additional processing is needed here
            vcf_content = [l.strip() for l in marker_hit_vcf_fh.readlines()
                           if re.match("##", l) is None]
            vcf_content_tab = [row_string.split("\t") for row_string in vcf_content]
            marker_hit_vcf_df = pd.DataFrame(vcf_content_tab[1:],
                                          columns=vcf_content_tab[0])
            marker_hit_vcf_df = marker_hit_vcf_df.rename(columns = {"#CHROM":"CHROM"})
            
        return marker_hit_vcf_df
    
    def query_gnomad(chrom, pos, ref, alt):
        # refer to https://discuss.gnomad.broadinstitute.org/t/blocked-when-using-api-to-get-af/149/2
        api_url = 'https://gnomad.broadinstitute.org/api'
        gnomad_v = "gnomad_r4"
    
        query = '''
        query ($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                genome {
                    af
                }
            }
        }
        '''
        variant_id = f'{chrom}-{pos}-{ref}-{alt}'
        variables = {"variantId": variant_id, "dataset": gnomad_v}
    
        response = requests.post(
            api_url,
            json={'query': query, 'variables': variables})
        return response.json()
    
    def identify_af(chrom, pos, ref, alt):
        max_try_num = 5
        try_ix = 0
        fetch_response = dict()
        while (try_ix < max_try_num) & (len(fetch_response) == 0):
            try:
                time.sleep(0.5)
                fetch_response = query_gnomad(chrom, pos, ref, alt)
            except:
                # Access to the API is rate-limited to 10 requests per IP address per 60-second period.
                print("Fetcthing failed ({})...retry in 60 sec".format(try_ix))
                time.sleep(60)
                try_ix = try_ix+1
            
        if try_ix == max_try_num:
            warnings.warn("variant not found")
            return -1
        elif "errors" in fetch_response:
            warnings.warn("variant not found")
            return -1
        else:
            try:
                af = fetch_response["data"]["variant"]["genome"]["af"]
                print("SUCCCESSFUL FETCHING")
                return af
            except TypeError as e:
                warnings.warn("successful fetching, but no relevant data found")
                print(fetch_response)
                print(e)
                return -1

    marker_hit_vcf_df = read_vcf(vcf_f)
    af_series = marker_hit_vcf_df.progress_apply(lambda row: identify_af(
        int(row["CHROM"]), int(row["POS"]), str(row["REF"]), str(row["ALT"])),
        axis=1)
    af_series = af_series.rename("AF")
    out_df = pd.concat([marker_hit_vcf_df[["CHROM", "POS", "REF", "ALT"]],
                        af_series], axis=1)
    return out_df

vcf_f = "scripts/annotation_clnvar_matched.vcf"
standard_out_df = standard_pipeline(vcf_f)
standard_out_df = standard_out_df.sort_values("AF")

out_dir = "scripts/test_af/"
standard_out_f = "standard_af_test.txt"
standard_out_f = os.path.join(out_dir, standard_out_f)
with open (standard_out_f, "w") as standard_out_fh:
    standard_out_df.to_csv(standard_out_fh, header=True, index=None, sep = "\t")

test_out_f = "test_af_test.txt"
test_out_f = os.path.join(out_dir, test_out_f)
test_out_df = GnomadUtils.make_af_table(vcf_f=vcf_f, out_f=test_out_f)
test_out_df = test_out_df.sort_values("AF")
assert os.path.isfile(test_out_f) & (os.path.getsize(test_out_f) > 0)
assert test_out_df.equals(standard_out_df)
print("DONE")   