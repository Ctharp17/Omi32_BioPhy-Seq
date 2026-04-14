# Omi32_BioPhy-Seq

**GitHub repository for analysis and figure generation in Tharp et al 2026 "Biophysical trade-offs in antibody evolution are resolved by conformation-mediated epistasis"**

All of the analysis was performed in python and most of the analysis uses the conda environment provided in `environment.yml` at the root of this directory (same working directory as this README file). This environment was exported with the `--from-history` flag, so it should work across Linux and MacOS systems. To install it, use `conda env create -f environment.yml` and activate it with `conda activate Tharp-et-al-2026`. 

Several analyses (i.e., barcode-variant-table generation, BioPhy-Seq high-throughput measurements, force-directed-layouts, and Figure 5/S-Figure 9) require separate conda environments, which are specified below and in READMEs in the requisite analysis folders.

The repository is separated into several folders:

**Figures**
- This folder contains the code used to generate Figures 2-5 (main text) and Supplemental Figures 1-8.
- All file paths for these scripts (input data) are directly linked to upstream analyses in the other folders described below.
- Some figures are generated in analysis folders for simplicity, for any figure that this is the case (e.g., Figure 3c) a README file is provided in that figure's subfolder (e.g, `Figures/Figure_3/README.md`) describing where those figure files can be found/generated.
- All code in `Figures` can be run with the global environment `Tharp-et-al-2026` except for code in `Figure_5-S_Figure_9`; a `README` file and environment file are provided in that subfolder to run the analysis.

**barcode-variant-table_generation**
- `link_genos_barcodes` contains the Bloom Lab's `dms-vep-pipeline-3`, which was used to analyze circular consensus sequences from PacBio long-read sequencing data and assign DNA-barcodes to Omi32 combinatorial library antibody variants. This subfolder comes with its own `README` and environment file for running the analysis, along with the necessary configuration and input files. We are still uploading the raw fastq files from the PacBio sequencing run and will update this section with the link to those data when available.
- `clean_geno-barcode_table` contains a custom script for cleaning the `dms-vep-pipeline-3` results such that we get a barcode-variant table with only the correct (designed) Omi32 mutations and DNA barcodes for the high-throughput BioPhy-Seq analysis (see below). The analysis in this folder can be run with the global environment `Tharp-et-al-2026`.

**BioPhy-Seq_measurements**
- `high-throughput_measurements` contains our BioPhy-Seq high-throughput measurement inference pipeline, which was used to infer equilibrium binding affinities to SARS-CoV-2 spike RBDs, antibody expression, and polyspecificity. This subfolder comes with its own `README` and environment file for running the analysis, along with the necessary configuration and input files. We are still uploading the raw fastq files from the Illumina sequencing runs and will update this section with the link to those data when available. See `results_Omi-32` in this subfolder for the raw measurement outputs from the pipeline and `results_Omi-32/cleaned_datasets` for the final 'cleaned' datasets (and the code used to filter the raw outputs).
- `isogenic_low-throughput_measurements` contains low-throughput validation measurements (shown in Supplemental Figure 1) inferred directly from flow cytometry data. All analyses in this folder can be run with the global environment `Tharp-et-al-2026`.

**epistasis_inference**
- `datasets_for_inference` contains a script for parsing the cleaned measurement files to assign mutation location positions and provide a final QC check of the data and binary genotype encodings (we strongly recommend NOT opening and saving files with the binary genotype (`geno`) column, as spreadsheet managers will often remove leading zeros or convert to scientific notation and remove lagging zeros which can cause major artifacts when inferring epistatic interactions between mutations)
- `global_epistasis_inference` contains subfolders and scripts for applying an initial nonlinear transformation to all of the BioPhy-Seq datasets for both reference frameworks (i.e., reference-based and reference-free) - tries out multiple nonlinear functions for fitting the data and prints out a comparison of fit qualities for each function, framework, and high-throughput dataset
- `linear_interaction_models` contains subfolders and scripts for epistasis inference on all of BioPhy-Seq datasets. This includes the reference-free and reference-based approaches, with and without initial non-linear transformation, cross-validation of model performance, and final outputs (i.e., model performance, epistatic coefficients, and heatmaps for each individual model tested). **The predicted phenotypes used in Figures 2-3 and Supplemental Figures 1-4 can be found under the reference-based subfolder for each measurement**
- All analysis in this folder can be run with the global environment `Tharp-et-al-2026`.

**pathway_analysis**
- This folder contains the pathway analysis performed for Figure 3 and Supplemental Figure 4. The notebooks that run this analysis are modular (i.e., most of the parameters in the configuration section of the notebook can be changed and easily rerun), but have been saved as separate notebook runs for simplicity.
- Each of the models (affinity-only, antigen-capture, and competitive-capture) have their own subfolders, and each subfolder contains the results from 'strong' and 'weak' selection runs
- All analyses in this folder can be run with the global environment `Tharp-et-al-2026`

**force-directed-layouts**
- This folder contains the code used to generate force-directed layouts of the Omi32 biophysical data and the interactive web browser for this analysis (published at https://amphilli.github.io/Omi32_browser_git/).
- While we only show the force-directed layout generated from the BA.1 affinity data and overlay the other datasets on top of it, this code can be used to generate layouts with weights calculated from the other datasets as well.
- This subfolder comes with its own `README` and environment file for running the analysis, along with the necessary configuration and input files.
