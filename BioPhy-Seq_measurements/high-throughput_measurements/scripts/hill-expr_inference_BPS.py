### infer expression, Kds, and/or PSR EC50 measurements with fit parameters set in config file
## can specify in config whether you want to do a three or four parameter fit (4 parameter fit allows you to try different hill coefficients if desired)
   ## we chose a 3 parameter fit for Kd (assuming a 1-site binding model) and a 3 parameter fit for PSR EC50 (best fit our data)

import pandas as pd
import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
mpl.rcParams['figure.dpi'] = 300

conf = snakemake.config
fit_hill = conf.get('fit_hill', False)  # True => 4p Hill; False => original 3p

# bring in initial (guess) parameters for hill fit from config file (different for PSR and Kd!)
def _p0_3p():
    d = conf.get('p0_3p', {})
    return [
        float(d.get('logEC50', -9.0)),
        float(d.get('A', 1e4)),
        float(d.get('B', 1e2)),
    ]

# bring in bound parameters for hill fit from config file (different for PSR and Kd!)
def _bounds_3p():
    b = conf.get('bounds_3p', {})
    lower = [
        float(b.get('logEC50_min', -12)),
        float(b.get('A_min', 1e2)),
        float(b.get('B_min', 1e1)),
    ]
    upper = [
        float(b.get('logEC50_max', -1)),
        float(b.get('A_max', 1e6)),
        float(b.get('B_max', 1e5)),
    ]
    return (lower, upper)

# bring in maxfev (number of fit iterations to try to minimize residuals)
def _maxfev_3p():
    return int(conf.get('maxfev_3p', 400000))

# same as above but for 4 parameter fit - only difference here is that there are additional guess and bounding parameters for the Hill Coefficient (N)
def _p0_4p():
    d = conf.get('p0_4p', {})
    return [
        float(d.get('logEC50', -9.0)),
        float(d.get('A', 1e4)),
        float(d.get('B', 1e2)),
        float(d.get('n', 1.0)),
    ]

def _bounds_4p():
    b = conf.get('bounds_4p', {})
    lower = [
        float(b.get('logEC50_min', -12)),
        float(b.get('A_min', 1e2)),
        float(b.get('B_min', 1e1)),
        float(b.get('n_min', 0.2)),
    ]
    upper = [
        float(b.get('logEC50_max', -1)),
        float(b.get('A_max', 1e6)),
        float(b.get('B_max', 1e5)),
        float(b.get('n_max', 5.0)),
    ]
    return (lower, upper)

def _maxfev_4p():
    return int(conf.get('maxfev_4p', 1000000))

# load sample info
sample_info = snakemake.params.sample_info
construct = sample_info.construct.unique()[0]
replicate = sample_info.replicate.unique()[0]

# use 'solute' identity to decide what to call the hill fit result in the columns
solute = None
if 'solute' in sample_info.columns and not sample_info['solute'].empty:
    solute = sample_info['solute'].iloc[0]
kd_col = "PSR_log10EC50" if (solute == "PSR") else "log10Kd"

df = pd.read_csv(snakemake.input.count_table, sep="\t")
df["geno"] = df.geno.apply(lambda x: f"{int(x):0{conf['num_mutations']}d}")

fluor = pd.read_csv(snakemake.input.fluorescence, sep="\t")
fluor = fluor[fluor.replicate == replicate]

# calculate mean expression (if data exists for this run - see methods for easier to read mathematical explanation) 
def mean_expression(sample_info, fluor, df):
    sample_info = sample_info[sample_info.concentration == 'F'].copy()
    fluor = fluor[fluor.concentration == 'F'].copy()
    nb_bins = sample_info.bin.nunique()
    nb_genos = len(df)
    probas = np.zeros((nb_bins, nb_genos))
    counts = np.zeros((nb_bins, nb_genos))
    cells = np.zeros(nb_bins)
    meanfluor, stdfluor = np.zeros((2, nb_bins))

    for bb, gate in enumerate(range(1, nb_bins+1)):
        col_name = f"{construct}_{replicate}_F_{gate}"
        if col_name not in df.columns:
            continue
        counts[bb, :] = df[col_name].fillna(0).to_numpy()

        cell_info = sample_info[sample_info.bin == gate]
        if not cell_info.empty:
            cells[bb] = cell_info["cell count"].iloc[0]

        fluor_info = fluor[fluor.bin == gate]
        if not fluor_info.empty:
            meanfluor[bb] = fluor_info["mean_log10_FITCs"].iloc[0]
            stdfluor[bb] = fluor_info["std_log10_FITCs"].iloc[0]

    with np.errstate(divide='ignore', invalid='ignore'):
        total_counts_per_bin = counts.sum(axis=1)
        valid_bins = (total_counts_per_bin > 0) & (cells > 0)
        safe_total = np.where(valid_bins, total_counts_per_bin, 1)
        scaled_counts = np.divide(counts, safe_total[:, None]) * cells[:, None] * valid_bins[:, None]

        probas_sum = scaled_counts.sum(axis=0)
        probas = np.divide(
            scaled_counts,
            probas_sum[None, :],
            out=np.zeros_like(scaled_counts),
            where=probas_sum[None, :] != 0
        )

        mean_log10_fluor = (probas * meanfluor[:, None]).sum(axis=0)
        std_log10_fluor = np.sqrt((stdfluor[:, None]**2 * probas**2 +
                                   np.divide(meanfluor[:, None]**2 * probas**2,
                                             counts,
                                             out=np.zeros_like(probas),
                                             where=counts != 0)
                                  ).sum(axis=0))
        # clip st_log_10 fluor zeros to 1e-4 so they aren't dropped later - had an issue where antibody variants that only fall into one bin end up with std_fluor of 0 resulting in divide by zero errors when fitting curves and calculating error
        std_log10_fluor = np.clip(std_log10_fluor, a_min=1e-4, a_max=None)

    return mean_log10_fluor, std_log10_fluor

# 3-parameter hill function (one-site binding model) used for all inferences in this study
def sigmoid_3param(c, logEC50, A, B):
    """
    3-parameter model (n = 1), same shape as your original:
        log10F = log10( A * (10^c / (10^c + 10^logEC50)) + B )
    """
    return np.log10(A * (10.0**c / ((10.0**c) + (10.0**logEC50))) + B)

# 4-parameter hill function (can infer cooperativity > or < 1)
def hill_full_logF(c, logEC50, A, B, n):
    """
    Expanded, easy-to-read 4-parameter Hill model (no algebraic shortcuts):
        F_linear = A * ((10^c)^n / ((10^c)^n + (10^logEC50)^n)) + B
        return log10(F_linear)
    """
    L_linear    = 10.0 ** c
    EC50_linear = 10.0 ** logEC50
    L_pow    = L_linear ** n
    EC50_pow = EC50_linear ** n
    frac = L_pow / (L_pow + EC50_pow)
    F_linear = A * frac + B
    return np.log10(F_linear)


def extractKd_3param(concentrations, bins_mean, bins_std):
    p0 = _p0_3p()
    bounds = _bounds_3p()
    popt, pcov = scipy.optimize.curve_fit(
        sigmoid_3param,
        concentrations, bins_mean,
        p0=p0,
        sigma=bins_std, absolute_sigma=True,
        bounds=bounds,
        maxfev=_maxfev_3p()
    )
    # Preserve sign convention: store +log10Kd (or PSR_log10EC50) in df
    log10Kd_est = -1 * popt[0]
    yfit = sigmoid_3param(concentrations, *popt)
    r2 = 1 - np.sum((yfit - bins_mean)**2) / np.sum((bins_mean - bins_mean.mean())**2)
    sigma_log10Kd = np.sqrt(pcov[0, 0]) if (pcov.shape == (3, 3) and pcov[0, 0] >= 0) else np.nan
    return (log10Kd_est, popt[1], popt[2], r2, sigma_log10Kd)

def extractKd_4param(concentrations, bins_mean, bins_std):
    p0 = _p0_4p()
    bounds = _bounds_4p()
    popt, pcov = scipy.optimize.curve_fit(
        hill_full_logF,
        concentrations, bins_mean,
        p0=p0,
        sigma=bins_std, absolute_sigma=True,
        bounds=bounds,
        maxfev=_maxfev_4p()
    )
    # Preserve same convention for the primary column
    log10Kd_est = -1 * popt[0]
    yfit = hill_full_logF(concentrations, *popt)
    r2 = 1 - np.sum((yfit - bins_mean)**2) / np.sum((bins_mean - bins_mean.mean())**2)
    sigma_log10Kd = np.sqrt(pcov[0, 0]) if (pcov.shape == (4, 4) and pcov[0, 0] >= 0) else np.nan

    hill_n = float(popt[3])
    hill_n_sigma = np.sqrt(pcov[3, 3]) if (pcov.shape == (4, 4) and pcov[3, 3] >= 0) else np.nan
    return (log10Kd_est, popt[1], popt[2], r2, sigma_log10Kd, hill_n, hill_n_sigma)

#compute Kds/PSR EC50s for datasets, saving fit parameters (with error) along with median sequencing coverage for each antibody variant across each concentration
def compute_Kds(sample_info, fluor, df_all):
    sample_info = sample_info[~sample_info.concentration_float.isna()].copy()
    nb_bins = sample_info.bin.nunique()
    nb_concs = sample_info.concentration.nunique()
    concentrations = sample_info.concentration_float.unique()
    nb_genos = len(df_all)
    probas = np.zeros((nb_bins, nb_genos, nb_concs))
    counts = np.zeros((nb_bins, nb_genos, nb_concs))
    cells = np.zeros((nb_bins, nb_concs))
    meanfluor, stdfluor = np.zeros((2, nb_bins, nb_concs))
    for bb, gate in enumerate(range(1, nb_bins+1)):
        for cc, conc in enumerate(sample_info.concentration.unique()):
            col_name = f"{construct}_{replicate}_{conc}_{gate}"
            if col_name not in df_all.columns:
                continue
            counts[bb, :, cc] = df_all[col_name].fillna(0).to_numpy()
            cell_info = sample_info[(sample_info.concentration == conc) & (sample_info.bin == gate)]
            if not cell_info.empty:
                cells[bb, cc] = cell_info["cell count"].iloc[0]
            fluor_info = fluor[(fluor.concentration == conc) & (fluor.bin == gate)]
            if not fluor_info.empty:
                meanfluor[bb, cc] = fluor_info["mean_log10_PEs"].iloc[0]
                stdfluor[bb, cc] = fluor_info["std_log10_PEs"].iloc[0]

    with np.errstate(divide='ignore', invalid='ignore'):
        total_counts_per_bin = counts.sum(axis=1)
        valid_bins = (total_counts_per_bin > 0) & (cells > 0)
        safe_total = np.where(valid_bins, total_counts_per_bin, 1)
        scaled_counts = np.divide(counts, safe_total[:, None, :]) * cells[:, None, :] * valid_bins[:, None, :]

        proba_sum = scaled_counts.sum(axis=0)
        probas = np.divide(
            scaled_counts, proba_sum[None, :, :],
            out=np.zeros_like(scaled_counts),
            where=proba_sum[None, :, :] != 0
        )

        mean_log10_fluor = (probas * meanfluor[:, None, :]).sum(axis=0)
        std_log10_fluor = np.sqrt((stdfluor[:, None, :]**2 * probas**2 +
                                   np.divide(meanfluor[:, None, :]**2 * probas**2,
                                             counts,
                                             out=np.zeros_like(probas),
                                             where=counts != 0)
                                  ).sum(axis=0))
        # Clip zeros to 1e-4 (see explanation above)
        std_log10_fluor = np.clip(std_log10_fluor, a_min=1e-4, a_max=None)

        total_cells_valid_per_conc = (cells * valid_bins).sum(axis=0)          
        variant_cells_per_conc = scaled_counts.sum(axis=0)                     
        variant_frac_per_conc = np.divide(
            variant_cells_per_conc,
            total_cells_valid_per_conc[None, :],
            out=np.full_like(variant_cells_per_conc, np.nan, dtype=float),
            where=total_cells_valid_per_conc[None, :] != 0
        )
        median_variant_fraction = np.nanmedian(variant_frac_per_conc, axis=1)  

    # expose via module-level globals
    global _variant_frac_per_conc, _median_variant_fraction
    _variant_frac_per_conc = variant_frac_per_conc
    _median_variant_fraction = median_variant_fraction

    concentrations = -concentrations

    if fit_hill:
        Kds, A, B, err, cov, hill_n, hill_n_sigma = np.zeros((7, nb_genos))
    else:
        Kds, A, B, err, cov = np.zeros((5, nb_genos))

    for s in range(nb_genos):
        notnanindex = [ii for ii in range(nb_concs)
                       if not np.isnan(mean_log10_fluor[s, ii] + std_log10_fluor[s, ii]) and mean_log10_fluor[s, ii] != 0]
        if len(notnanindex) < 4:
            if fit_hill:
                Kds[s], A[s], B[s], err[s], cov[s], hill_n[s], hill_n_sigma[s] = [np.nan]*7
            else:
                Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        elif np.sum(counts.sum(axis=0) > conf['min_number_counts']) < 4:
            if fit_hill:
                Kds[s], A[s], B[s], err[s], cov[s], hill_n[s], hill_n_sigma[s] = [np.nan]*7
            else:
                Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        else:
            if fit_hill:
                (Kds[s], A[s], B[s], err[s], cov[s],
                 hill_n[s], hill_n_sigma[s]) = extractKd_4param(
                    concentrations[notnanindex],
                    mean_log10_fluor[s, notnanindex],
                    std_log10_fluor[s, notnanindex]
                )
            else:
                Kds[s], A[s], B[s], err[s], cov[s] = extractKd_3param(
                    concentrations[notnanindex],
                    mean_log10_fluor[s, notnanindex],
                    std_log10_fluor[s, notnanindex]
                )

    if fit_hill:
        return Kds, A, B, err, cov, hill_n, hill_n_sigma, mean_log10_fluor, std_log10_fluor, concentrations
    else:
        return Kds, A, B, err, cov, mean_log10_fluor, std_log10_fluor, concentrations


#plotting function for quick QC of a random and small subset of curve fits
xs = np.linspace(-14, -6, 100)
df_subset = df.sample(n=8)

if fit_hill:
    Kds, As, Bs, errs, covs, hill_ns, hill_ns_sigmas, mlog10, slog10, concentrations = compute_Kds(sample_info, fluor, df_subset)
else:
    Kds, As, Bs, errs, covs, mlog10, slog10, concentrations = compute_Kds(sample_info, fluor, df_subset)

xlim = (-max(sample_info.concentration_float)-0.5,
        -min(sample_info[sample_info.concentration_float != 0].concentration_float)+0.5)
fig, ax = plt.subplots(4, 2, figsize=(10, 10), sharex=True)
for ii, ax in enumerate(ax.flatten()):
    ax.errorbar(x=concentrations, y=mlog10[ii, :], yerr=slog10[ii, :], label=df_subset.geno.iloc[ii])
    if fit_hill:
        ax.plot(xs, hill_full_logF(xs, -Kds[ii], As[ii], Bs[ii], hill_ns[ii]) if np.isfinite(hill_ns[ii]) else np.nan)
    else:
        ax.plot(xs, sigmoid_3param(xs, -Kds[ii], As[ii], Bs[ii]))
    ax.set_xlim(xlim)
    ax.set_xlabel("$\\log_{10}(\\mathrm{Concentration})$")
    ax.set_ylabel("Est. Mean Fluorescence")
plt.savefig(snakemake.output.plot_test_curve)

#now infer all the Kds/PSR EC50s and save the data
if fit_hill:
    Kds, As, Bs, errs, covs, hill_ns, hill_ns_sigmas, mean_log10_PE, std_log10_PE, concs = compute_Kds(sample_info, fluor, df)
else:
    Kds, As, Bs, errs, covs, mean_log10_PE, std_log10_PE, concs = compute_Kds(sample_info, fluor, df)

df[kd_col] = Kds
df["A"] = As
df["B"] = Bs
df["r2"] = errs
df["sigma"] = covs

if fit_hill:
    df["hill_n"] = hill_ns
    df["hill_n_sigma"] = hill_ns_sigmas

for cc in range(mean_log10_PE.shape[1]):
    df[f"mean_log10PE{cc}"] = mean_log10_PE[:, cc]
    df[f"std_log10PE{cc}"]  = std_log10_PE[:, cc]  

df["median_variant_fraction"] = _median_variant_fraction
for cc in range(mean_log10_PE.shape[1]):
    df[f"variant_frac_conc{cc}"] = _variant_frac_per_conc[:, cc]

if 'F' in sample_info.concentration.unique():
    m10s, s10s = mean_expression(sample_info, fluor, df)
    df["Mean fluorescence expression"] = m10s
    df["Std fluorescence expression"]  = s10s

base_cols = ["geno", kd_col, "A", "B", "r2", "sigma", "median_variant_fraction"]
pe_cols   = [f"mean_log10PE{cc}" for cc in range(mean_log10_PE.shape[1])] + \
            [f"std_log10PE{cc}"  for cc in range(mean_log10_PE.shape[1])]
frac_cols = [f"variant_frac_conc{cc}" for cc in range(mean_log10_PE.shape[1])]

out_cols = base_cols + pe_cols + frac_cols
if 'F' in sample_info.concentration.unique():
    out_cols += ["Mean fluorescence expression", "Std fluorescence expression"]
if fit_hill:
    base_idx = out_cols.index("sigma") + 1
    out_cols = out_cols[:base_idx] + ["hill_n", "hill_n_sigma"] + out_cols[base_idx:]

df[out_cols].to_csv(snakemake.output.tsv, sep="\t", index=False)

# optional plots for assessing inference quality for each dataset - these can throw bugs so default is set to False
with pd.option_context('mode.use_inf_as_na', True):
    fig, ax = plt.subplots()
    if conf.get('do_plots', True):
        sns.histplot(x=kd_col, data=df.dropna(subset=[kd_col]),
                     ax=ax, label="All intermediates")
        if df.geno.apply(lambda x: '1' not in x).sum() > 0 and kd_col in df:
            ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)][kd_col].iloc[0], label="Original", c="g")
        if df.geno.apply(lambda x: '0' not in x).sum() > 0 and kd_col in df:
            ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)][kd_col].iloc[0], label="Variant", c="r")
        ax.legend()
    plt.savefig(snakemake.output.plot_Kd_distribution)

    fig, ax = plt.subplots()
    if "Mean fluorescence expression" in df and conf.get('do_plots', True):
        sns.histplot(x="Mean fluorescence expression", data=df.dropna(subset=["Mean fluorescence expression"]),
                     label="All intermediates", ax=ax)
        if df.geno.apply(lambda x: '1' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)]["Mean fluorescence expression"].iloc[0],
                       label="Original", c="g")
        if df.geno.apply(lambda x: '0' not in x).sum() > 0:
            ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)]["Mean fluorescence expression"].iloc[0],
                       label="Variant", c="r")
        ax.legend()
    plt.savefig(snakemake.output.plot_fluo_distribution)

    if "Mean fluorescence expression" in df and conf.get('do_plots', True):
        sns.jointplot(x="Mean fluorescence expression", y=kd_col,
                      data=df.dropna(subset=["Mean fluorescence expression", kd_col]),
                      kind="hex", bins='log')
    else:
        plt.subplots()
    plt.savefig(snakemake.output.plot_corr_fluo_Kd)

    if conf.get('do_plots', True):
        df["$\\log_{10}(\\mathrm{error})$"] = np.log(1 - df["r2"])
        sns.jointplot(x="$\\log_{10}(\\mathrm{error})$", y=kd_col,
                      data=df.dropna(subset=["$\\log_{10}(\\mathrm{error})$", kd_col]),
                      kind="hex", bins='log')
    else:
        plt.subplots()
    plt.savefig(snakemake.output.plot_corr_Kd_error)
