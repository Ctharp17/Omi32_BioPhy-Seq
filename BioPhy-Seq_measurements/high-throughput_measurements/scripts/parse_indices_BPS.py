### use regexes generated earlier to parse in-line indices and ensure each read from each sample (bin) has the correct set of in-line indices

import regex
import pandas as pd
import numpy as np

# compile fuzzy regexes for in-line indices based on allowed number of substitutions (from config file)
n_subs = snakemake.config['num_substitutions_tolerated_in_inline_idx']
col_idx_regex = regex.compile('(%s){s<=%d}' % (snakemake.params.col_inline_idx, n_subs))
row_idx_regex = regex.compile('(%s){s<=%d}' % (snakemake.params.row_inline_idx, n_subs))


nb_reads = 0
reads_thrown_out = 0
reads_with_no_index_err = 0
reads_with_corrected_index_err = 0

# align in-line indices --> if read 1 and read 2 indices are correct, save UMI and barcode information for downstream deduplication and counting
with open(snakemake.input.tsv) as f, open(snakemake.output.tsv, 'w') as fw:
    f.readline() 
    fw.write("\t".join(["UMI", "barcode"]) + "\n")
    for line in f:
        col_idx, row_idx, UMI, barcode = line.strip().split("\t")
        nb_reads += 1
        match_1 = col_idx_regex.match(col_idx)
        match_2 = row_idx_regex.match(row_idx)

        if match_1 is None or match_2 is None:
            reads_thrown_out += 1
            continue

        if match_1.fuzzy_counts[0] == 0 and match_2.fuzzy_counts[0] == 0:
            reads_with_no_index_err += 1
        else:
            reads_with_corrected_index_err += 1

        fw.write("\t".join([UMI, barcode]) + "\n")

# Write stats
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'reads_with_corrected_index_err' : [reads_with_corrected_index_err],
    'reads_with_no_index_err' : [reads_with_no_index_err],
    'reads_thrown_out' : [reads_thrown_out],
    'fraction_of_index_accepted' : [(reads_with_corrected_index_err +
                                     reads_with_no_index_err) / nb_reads
                                    ] if nb_reads != 0 else np.nan,
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)
