### generating regexes to enable parsing of Illumina paired-end reads
## regexes are used to check in-line index usage, align to each read's assigned constant region, and identify UMIs and barcodes

import math
from Bio.Seq import Seq
import pandas as pd
import regex
import numpy as np


def rc(x):
    return str(Seq(x).reverse_complement())

UMI_length = snakemake.config['umi_length']
col_inline_idx = snakemake.config['col_inline_indices']
row_inline_idx = snakemake.config['row_inline_indices']

num_idx_subs = snakemake.config['num_substitutions_tolerated_in_inline_idx']
num_substitutions_tolerated_per_ten_bp = snakemake.config['num_substitutions_tolerated_per_ten_bp']
are_indels_tolerated_in_primers = snakemake.config['are_indels_tolerated_in_primers']

#return an OR regex that can handle all the different in-line index combinations that we used for the multiplexed Illumina sequencing libraries
def make_or_regex(seqs, name, dist=0):
    ret = f"(?P<{name}>{'|'.join(seqs)})"
    if dist != 0:
        ret += '{s<=%d}' % dist
    return ret

#return fuzzy_regex that can handle set amounts of error in the read sequences (based off config file)
def make_fuzzy_regex(search, num_err, name='', bestmatch=True, allow_indels=False):
    if search == '':
        return ''
    r = ''
    if bestmatch and num_err != 0:
        r += '(?e)'
    r += '('
    if name != '':
        r += f'?P<{name}>'
    r += search
    r += ')'
    if num_err != 0:
        r += '{%c<=%d}' % ('e' if allow_indels else 's', num_err)
    return r

#assemble the regex components
def assemble_regex(const_region, num):
    regex_parts = []
    num_mut_tolerated = math.ceil(
        (len(const_region) / 10) * num_substitutions_tolerated_per_ten_bp)
    is_indel_tolerated = are_indels_tolerated_in_primers

    regex_parts.append(make_fuzzy_regex(const_region,
                                        num_err=num_mut_tolerated,
                                        allow_indels=is_indel_tolerated,
                                        name=f'const_{num}'))
    print('\n'.join(regex_parts))
    return ''.join(regex_parts)



read_1_const_1 = snakemake.config['for_const_reg_1']
read_1_const_2 = snakemake.config['for_const_reg_2']
read_2_const_1 = snakemake.config['rev_const_reg_1']

read_2_const_1_rc = rc(read_2_const_1)

#make regexes for all amplicon components
UMI_regex = f"(?P<UMI>[ACGT]{{{UMI_length}}})"
cols_regex = make_or_regex(col_inline_idx, name='index', dist=num_idx_subs)
rows_regex = make_or_regex(row_inline_idx, name='index', dist=num_idx_subs)
read_1_const_1_regex = assemble_regex(read_1_const_1,1)
read_1_const_2_regex = assemble_regex(read_1_const_2,2)
read_2_const_1_regex = assemble_regex(read_2_const_1_rc,3)
barcode_regex = "(?P<barcode>[ACGT]{2,30})"

#concatenate regexes for fastq parsing and UMI/barcode assignment
read_1_regex = cols_regex + read_1_const_1_regex + barcode_regex + read_1_const_2_regex
read_2_regex = UMI_regex + rows_regex + read_2_const_1_regex

with open(snakemake.output.read_1_regex, 'w') as f:
    f.write(read_1_regex)

with open(snakemake.output.read_2_regex, 'w') as f:
    f.write(read_2_regex)
