
# BioPhy-Seq analysis pipeline for inferring affinity to Wuhan-Hu, BA.1, and BA.4 SARS-CoV-2 spike RBD, antibody expression, and polyspecificity for the Omi-32 combinatorial library in Tharp et al, 2026.

The code here encompasses a Snakemake pipeline that can be used to run the analysis (read the docs here: https://snakemake.readthedocs.io/en/stable/).

This pipeline was initially developed in Phillips and Maurer et al, 2022 by Jeffrey Chang and was adapted for this study to incorporate barcoded antibody variants (see methods) and measure polyspecificity in addition to antibody surface expression and antigen affinity by Cole Tharp.

## About the pipeline

This analysis is done with the `snakemake` pipeline described in the [Snakefile](Snakefile). To run this, first make sure you have installed [conda](https://docs.conda.io/en/latest/) 


```

build the environment with `conda env create -f environment.yml` , please note that this environment configuration file is saved using the --from-history flag so it should enable automatic package handling between linux and MacOS.


```

Then activate conda with `conda activate BioPhy-Seq`


```

You can run the analysis with 'snakemake -j N' , Where you can replace `N` by the number of jobs you want to run at the same time (where 1 job usually utilizes 1 thread).


## How to run the analysis

### Configure the pipeline

1. Go over the settings in the files [metadata/config.yaml](metadata/config.yaml) about how to parse and analyze the data.
2. Enter information about each sample into [metadata/sample_info.tsv](metadata/sample_info.tsv).

### Prepare the data input

1. Enter the cell counts into [metadata/sample_info.tsv](metadata/sample_info.tsv).
2. Transfer the NGS read data to the correct input directory
    1. Place the files under `{data_dir}/fastq/`. If you performed separate sequencing runs, you can place each one in its own `data_dir`, as long as you specify it in the column of [metadata/sample_info.tsv](metadata/sample_info.tsv).
    2. Make sure the files are named `{sample_name}_S{sample_id}_L00?_R?_001.fastq.gz`, where `sample_name` and `sample_id` are the columns in [metadata/sample_info.tsv](metadata/sample_info.tsv).
3. Transfer the fcs files (flow cytometry information) to the correct input directory. Use Flow-Jo to read the fcs files. There's already a compensation applied. Select the "Sorted" tubes and export them to csv format (Export or Concatenate → CSV - Scale values, All compensated parameters). (You can also rerun the compensation with flow-Jo but that shouldn't be necessary.)
    1. Save the resulting files in `{data_dir}/{gate_events}`.
    2. Filenames should be in the format `export_Sorted_{construct}{replicate}_{concentration}_PE{bin}` for samples testing for binding and `export_Sorted_{construct}{replicate}_{concentration}_mycFITC_FITC{bin}` for the expression samples.

### Run snakemake

Snakemake will automatically perform the following steps:

1. Parse the fastq files into a table of the number of reads per genotype, for each the samples.
    1. Parse fastq files and separate each read into fields `col_idx, row_idx, UMI, barcode`
    2. Discard reads with incorrect inline indices, leaving fields `UMI, barcode`
    3. Use the variant-barcode dictionary (generated from PacBio long-read sequencing of antibody library) to assign binary genotypes to each barcode leaving fields 'UMI, barcode, geno'
    4. Count the number of unique UMIs appearing for each genotype, leaving fields `geno, count`
    5. Collect all samples into a single count table, with a column for each sample and a row for each genotype
2. Fit the Kd or PSR EC50 or mean fluorescence for each of the antibody variants
    1. Read the csv fluorescence files, correct for fluorescence below the linear range of the emission detector (on flow cytometer) and extract the mean and standard variance of the log-fluorescence.
    2. Calculate the mean bin from the cell counts, reads, and flourescence values
    3. Fit the one-site binding curve (can also specify to fit hill coefficient but default in a 3-parameter single-site binding model) to the mean bin.

The relevant output files will be the following:

* The count table is written to [results_{experiment_name}/count_table.tsv](results/count_table.tsv).
* The inferred Kds are written to `results_{experiment_name}/Kds/Kds_{construct}{replicate}.tsv'`
* The fluorescence information is written to [results_{experiment_name}/fluorescence.tsv](results/fluorescence.tsv).
* Statistics about sequencing errors, number of duplicate UMI's, etc. are written to the directory [results/stats/](results/stats)

**notes:

We did not run this analysis on a computing cluster, which is why no batch script for this analysis is included.

We ran this entire analysis pipeline on a custom-built linux (Ubuntu) machine with a AMD Ryzen Threadripper PRO 5975WX (32 cores, 64 threads), 2 TB SSD, and 128 GB RAM by running 'Snakemake -j30,' which took (on average) ~2 hours to run the entire analysis for biological replicates of psr or affinity measurements (if all of the fastq files are already stored and grouped locally). 


