# Pipeline adapted from Bloom Lab deep mutational scanning (DMS) of viral entry proteins (VEPs), version 3, to assign barcodes to antibody variants in the Omi32 combinatorial library from Tharp et al, 2026 (see README in dms-vep-pipeline-3 folder for more details on the full pipeline).

## Running the pipeline

1) Set working directory to the 'link_genos_barcodes' folder

2) install the conda environment from 'environment.yml' in this folder (this environment file should work for both macOS and linux platforms) with `conda env create -f environment.yml`

3) use `conda activate dms-vep-pipeline-3` to activate the environment

4) run the pipeline with `snakemake -jN -s dms-vep-pipeline-3/Snakefile`, where 'N' signifies the number of jobs you want to run at once (for this combinatorial library, we run N = 1 job at a time because the analysis performed by this pipeline for our specific usage isn't parallelizable)

## Important outputs and downstream usage of results

The output we care about for our study is 'results/variants/codon_variants_Omi-32.csv,' which contains each consensus barcode-to-antibody-variant relationship, along with read coverage stats. We further process this data to remove any non-designed mutations (from the Omi-32 combinatorial library) and any barcodes that have incorrect lengths. We then assign binary genotypes to each antibody variant, which are strings of 13 0/1s (where '0' denotes a wild-type/germline amino acid residue and '1' denotes a mutated/Omi32 amino acid residue). Code for this analysis is upstream of this folder at `../clean_geno-barcode_table/build_variant_table_with_binary_genos.ipynb`. The final output from this analysis is `../clean_geno-barcode_table/variants_barcodes_binary-genos.csv`. 

The `results` folder also contains notebooks (auto-generated from the Bloom Lab pipeline) that show statistics at each step of the process, from reading in and aligning PacBio circular consensus sequences (CCSs) to calling barcode-antibody-variant linkages.  