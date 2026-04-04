### calculate geometric mean fluorescence for each bin at each concentration for Kds/PSR/expression
## will need to change the column header in list_FITCs and list_PEs for whichever fluor/filter-laser configuration is used

import pandas as pd
import numpy as np

sample_table = snakemake.params.sample_table.copy()
sample_table = sample_table[sample_table.concentration != "unsorted"]

#change this based on your expression fluor/laser/filter configuration
def list_FITCs(csv):
    df = pd.read_csv(csv)
    try:
        return df["530_30 Blue B-A"].to_list()
    except KeyError:
        return df["FITC-A"].to_list()

#change this based on your Kd/PSR fluor/laser/filter configuration
def list_PEs(csv):
    df = pd.read_csv(csv)
    for col in ["582_15 YG D-A","PE_A"]:
        if col in df.columns:
            return df[col].to_list()
    return []

#pin raw values (not log-transformed yet) to a safe detection threshold within linear detector range (based on average CST settings provided by UCSF Flow Core)
def pin_below_threshold(values, threshold=100):
    return [max(v, threshold) for v in values]

# load the data
sample_table["FITCs"] = sample_table.flow_file.apply(list_FITCs)
sample_table["PEs"] = sample_table.flow_file.apply(list_PEs)

# pin values below linear range of cytometer emission detector
sample_table["pinned_FITCs"] = sample_table["FITCs"].apply(pin_below_threshold)
sample_table["pinned_PEs"] = sample_table["PEs"].apply(pin_below_threshold)

# log transform all fluorescence values
sample_table["log10FITCs"] = sample_table["pinned_FITCs"].apply(lambda vals: np.log10(np.array(vals)))
sample_table["log10PEs"] = sample_table["pinned_PEs"].apply(lambda vals: np.log10(np.array(vals)))

# compute the (geometric) mean and standard deviation
sample_table["mean_log10_PEs"] = sample_table.log10PEs.apply(np.mean)
sample_table["std_log10_PEs"] = sample_table.log10PEs.apply(np.std)
sample_table["mean_log10_FITCs"] = sample_table.log10FITCs.apply(np.mean)
sample_table["std_log10_FITCs"] = sample_table.log10FITCs.apply(np.std)

sample_table[["sample_id", "construct", "replicate", "concentration",
              "bin", "mean_log10_PEs", "std_log10_PEs",
              "mean_log10_FITCs", "std_log10_FITCs"]].to_csv(snakemake.output.tsv, sep="\t")
