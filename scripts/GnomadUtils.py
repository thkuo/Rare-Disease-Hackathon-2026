import requests
import time
import warnings
import pandas as pd
from tqdm import tqdm
from scripts import FileUtils

tqdm.pandas()

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

def make_af_table(vcf_f: str, out_f: str=None):
    vcf_df = FileUtils.read_vcf(vcf_f)
    af_series = vcf_df.progress_apply(lambda row: identify_af(
        int(row["CHROM"]), int(row["POS"]), str(row["REF"]), str(row["ALT"])),
        axis=1)
    af_series = af_series.rename("AF")
    out_df = pd.concat([vcf_df[["CHROM", "POS", "REF", "ALT"]],
                        af_series], axis=1)
    if out_f is not None:
        print("Saving result...")
        with open (out_f, "w") as out_fh:
            out_df.to_csv(out_fh, header=True, index=None, sep = "\t")

    return out_df