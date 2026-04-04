### assign genotypes to barcodes from Illumina reads

import regex
import pandas as pd
from Bio.Seq import Seq
import numpy as np
import gzip
import pybktree
import Levenshtein

# build a barcode:binary genotype dictionary 
variants_barcodes = pd.read_csv(snakemake.config["variant_barcode_table"], dtype = 'str')
variants_barcodes['barcode'] = variants_barcodes['barcode'].str.upper()

geno_binary = dict(zip(variants_barcodes['barcode'], variants_barcodes['binary_genotypes']))

# BK-tree builder for efficient searching of mismatched barcodes (https://github.com/benhoyt/pybktree)
## we do NOT implement BK-trees for final analyses (see config file) but did attempt to use them to save reads and improve data quality. For our data, specifically, allowing errors in barcodes did not improve our replicate correlations (and in some cases, made them worse)
barcodes = list(geno_binary.keys())
max_bc_err = snakemake.config['num_substitutions_tolerated_in_barcode']
bc_indels_tol = snakemake.config['allow_indels_in_barcode']

bk_tree = None
if max_bc_err > 0:
    ##BK-tree distance function
    def dist(a, b):
        return Levenshtein.distance(a, b)
    #build a BK-tree for all barcode entries
    bk_tree = pybktree.BKTree(dist, barcodes)

class BarcodeDoesntMatchError(ValueError):
    pass

#BK-tree lookup function
def mutations_to_binary_geno(barc):
    barc = str(barc).upper()

    # Exact match lookup
    if barc in geno_binary:
        return geno_binary[barc], barc, (0, 0, 0)

    # Skip fuzzy match if errors aren't allowed in barcode
    if max_bc_err == 0:
        raise BarcodeDoesntMatchError

    # BK-tree lookup
    candidates = bk_tree.find(barc, max_bc_err)
    if not candidates:
        raise BarcodeDoesntMatchError

    # Select the nearest neighbor
    best_dist, best_bc = min(candidates, key=lambda x: x[0])

    # Compute detailed edit operations
    ops = Levenshtein.editops(barc, best_bc)
    subs = sum(1 for op in ops if op[0] == 'replace')
    ins  = sum(1 for op in ops if op[0] == 'insert')
    dels = sum(1 for op in ops if op[0] == 'delete')

    return geno_binary[best_bc], best_bc, (subs, ins, dels)

incorrect_barcode = 0
reads_with_no_err = 0
reads_with_one_sub = 0
reads_with_two_plus_subs = 0
reads_with_indels = 0
total_accepted = 0
nb_reads = 0

with \
     open(snakemake.input.tsv, 'r') as f, \
     open(snakemake.output.tsv, 'w') as fw, \
     gzip.open(snakemake.output.bad_reads, 'wt') as fbad:

    f.readline() 
    fw.write("UMI\tgeno\tbarcode\n")
    fbad.write("barcode\n")

    for line in f:
        UMI, barcode = line.strip().split("\t")
        nb_reads += 1

        try:
            geno, bc, (subs_bc, ins_bc, dels_bc) = mutations_to_binary_geno(barcode)
        except BarcodeDoesntMatchError:
            incorrect_barcode += 1
            fbad.write(f"{barcode}\n")
            continue
        if not bc_indels_tol and (ins_bc + dels_bc) > 0:
            reads_with_indels += 1
            fbad.write(f"{barcode}\n")
            continue

        if subs_bc == 0 and (ins_bc + dels_bc) == 0:
            reads_with_no_err += 1
        elif (ins_bc + dels_bc) > 0:
            reads_with_indels += 1
        elif subs_bc == 1:
            reads_with_one_sub += 1
        else:
            reads_with_two_plus_subs += 1
        


        total_accepted += 1
        fw.write(f"{UMI}\t{geno}\t{bc}\n")


stats = pd.DataFrame({
    'sample': [snakemake.wildcards.sample],
    'reads_with_incorrect_barcode': [incorrect_barcode],
    'reads_with_no_err_in_barcode': [reads_with_no_err],
    'reads_with_one_sub_in_barcode': [reads_with_one_sub],
    'reads_with_two_or_more_subs_in_barcode': [reads_with_two_plus_subs],
    'reads_with_indels_in_barcode': [reads_with_indels],
    'total_accepted': [total_accepted],
    'fraction_accepted': [total_accepted/nb_reads if nb_reads else np.nan],
}).set_index('sample')

stats.to_csv(snakemake.output.stats, sep='\t')