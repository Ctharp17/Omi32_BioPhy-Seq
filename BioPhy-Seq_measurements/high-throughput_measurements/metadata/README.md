# Settings and sequences

These files are for processing sequencing and flow cytometry data for inference of Kds, PSR EC50s, and antibody expression in Tharp et al, 2026.


`sample_info.tsv`
-----------------
- `concentration`: Concentration of the solute, string for tube number/label
- `concentration_float`: Measured concentration of the solute in units of -log10 (Molar)

SPECIAL NOTE: since we can't take the log of zero (and for better accounting of background noise in flow cytometry measurements for our curve fits) we set our 'zero' (no antigen) concentration/sample to a number a couple of logs lower than the minimum concentration used (for this study, we used 10**(-14) Molar)

- `bin`: bin name used to label the tubes

**notes

Each set of biological duplicates has it's own sample_sheet.tsv and replicate_info.tsv (make sure to specify these in the appropriate config file if rerunning this analysis)

There are two config files, designated for affinity measurements and PSR EC50 measurements. These config files are identical, except for the single-site binding model fit parameters (i.e., initial guess parameters and concentration bounds) Please specify in the Snakefile which config file you are using.

