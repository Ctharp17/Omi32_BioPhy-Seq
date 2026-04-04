#!/usr/bin/env python3
"""
Epistasis Inference + Prediction Script

For each row in CONFIG_CSV, fits epistasis models (orders 1..max_order) and:
  1. Exports coefficients, covariance matrix, and metadata to MODEL_COEFFS_DIR
  2. If the row's 'predict_phenotypes' column is TRUE, also exports full
     2^L genotype predictions with uncertainty to PREDICTIONS_DIR (used for Figures 1-2 and ED Figures 1-2 in Tharp et al, 2026)

Expected columns in CONFIG_CSV:
  - scope              : "full" or "pos4"
  - phenotype          : "raw" or e.g. "latent_hill", "latent_logistic"
  - order              : max epistatic order to fit (int); fits 1..order inclusive
  - predict_phenotypes : TRUE/FALSE — whether to export predictions for this model

Latent CSV filename format:
  {ANTIGEN}_epistasis_linearized_{link_function}_v{VERSION}_{ENCODING}.csv
  where VERSION = 1 for scope="full", VERSION = 2 for scope="pos4"
"""

import os
import csv
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

# Configuration
ANTIGEN       = "ba1"
ENCODING      = "stat"  # "biochem" => 0/1 (reference-based encoding) or "stat" => -1/+1 (reference-free encoding)

CONFIG_CSV    = f"{ANTIGEN}_ep_orders_{ENCODING}.csv"
EPISTASIS_CSV = f"{ANTIGEN}_epistasis.csv"

GENO_COL  = "geno"
PHENO_COL = f"{ANTIGEN}_log10Kd_mean"
POS4_COL  = "pos4"

# Directory configuration

_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT        = os.path.normpath(os.path.join(_SCRIPT_DIR, "../../.."))
_ENC_SUBDIR  = "reference-based" if ENCODING == "biochem" else "reference-free"

RAW_DATA_DIR     = os.path.join(_ROOT, "datasets_for_inference")
LATENT_DATA_DIR  = os.path.normpath(os.path.join(_SCRIPT_DIR, "../../../global_epistasis_inference", ANTIGEN, _ENC_SUBDIR))

MODEL_COEFFS_DIR = os.path.join(_SCRIPT_DIR, "model_coeffs")
PREDICTIONS_DIR  = os.path.join(_SCRIPT_DIR, "predicted_phenotypes")

ALPHA          = 0.05  # for Bonferroni-corrected confidence intervals
POS4_MUT_INDEX = 3     # zero-based index of pos4 mutation (S75G)


# Helpers: data loading
def load_and_normalize_genos(df: pd.DataFrame, geno_col: str):
    df = df.copy()
    df[geno_col] = df[geno_col].astype(str).str.strip().str.replace(r"[^01]", "", regex=True)
    df = df.loc[df[geno_col].str.len() > 0].copy()
    L = int(df[geno_col].str.len().max())
    df[geno_col] = df[geno_col].str.zfill(L)
    return df, L


def encode_genotype_series(geno_series: pd.Series, encoding: str) -> np.ndarray:
    X = np.array([[int(c) for c in g] for g in geno_series], dtype=int)
    if encoding == "biochem":
        return X.astype(float)
    elif encoding == "stat":
        return (2 * X - 1).astype(float)
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


def encode_genotype_single(geno_str: str, encoding: str) -> np.ndarray:
    x = np.array([int(c) for c in geno_str], dtype=float)
    if encoding == "biochem":
        return x
    elif encoding == "stat":
        return 2 * x - 1
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


# Helpers: model fitting
def fit_epistasis_model(X: np.ndarray, y: np.ndarray, order: int,
                        feature_labels: list) -> dict:
    """
    Fit OLS epistasis model with PolynomialFeatures (interaction_only, include_bias).
    Returns dict with model, feature names, R², n_params, and fitted poly transformer.
    """
    poly = PolynomialFeatures(degree=order, include_bias=True, interaction_only=True)
    X_poly = poly.fit_transform(X)
    model = sm.OLS(y, X_poly).fit()
    feature_names = poly.get_feature_names_out(input_features=feature_labels)
    return {
        "model":         model,
        "feature_names": feature_names,
        "r2":            float(model.rsquared),
        "n_params":      int(len(model.params)),
        "poly":          poly,
    }


# Helpers: export
def export_coefficients(result: dict, output_path: str, alpha: float = 0.05) -> None:
    model         = result["model"]
    feature_names = result["feature_names"]
    n_params      = result["n_params"]
    alpha_corrected = alpha / float(n_params)
    conf_int = model.conf_int(alpha=alpha_corrected)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Params:",      n_params])
        writer.writerow(["Performance_R2:", result["r2"]])
        writer.writerow(["Alpha:",       alpha])
        writer.writerow(["Alpha_Bonf:",  alpha_corrected])
        writer.writerow(["Term", "Coefficient", "StdErr", "p_value", "CI_lower", "CI_upper"])
        for i in range(len(model.params)):
            term     = "Intercept" if i == 0 else ",".join(str(feature_names[i]).split(" "))
            coef     = float(model.params[i])
            stderr   = float(model.bse[i])
            pval     = float(model.pvalues[i])
            ci_lower = float(conf_int[i][0])
            ci_upper = float(conf_int[i][1])
            writer.writerow([term, coef, stderr, pval, ci_lower, ci_upper])

    print(f"  Coefficients : {output_path}")
    print(f"  R²={result['r2']:.4f}, n_params={n_params}, alpha_bonf={alpha_corrected:.3g}")


def export_covariance(result: dict, output_path: str) -> None:
    cov_df = pd.DataFrame(
        result["model"].cov_params(),
        index=result["feature_names"],
        columns=result["feature_names"],
    )
    cov_df.to_csv(output_path)
    print(f"  Covariance   : {output_path}  {result['model'].cov_params().shape}")


def export_metadata(result: dict, output_path: str, scope: str, phenotype: str,
                    order: int, encoding: str, n_samples: int,
                    feature_mask=None) -> None:
    metadata = {
        "scope":         scope,
        "phenotype":     phenotype,
        "order":         order,
        "encoding":      encoding,
        "n_params":      result["n_params"],
        "n_samples":     n_samples,
        "r2":            result["r2"],
        "feature_names": result["feature_names"].tolist(),
    }
    if feature_mask is not None:
        metadata["dropped_features"] = [i for i, m in enumerate(feature_mask) if not m]
    with open(output_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"  Metadata     : {output_path}")


# Helpers: prediction
def generate_all_genotypes(length: int) -> list:
    return [format(i, f"0{length}b") for i in range(2 ** length)]


def predict_with_uncertainty(geno_str: str, beta: np.ndarray, cov_beta: np.ndarray,
                             poly: PolynomialFeatures, encoding: str,
                             dropped_features: list = None) -> tuple:
    x_base = encode_genotype_single(geno_str, encoding)
    if dropped_features:
        mask = np.ones(len(x_base), dtype=bool)
        mask[dropped_features] = False
        x_base = x_base[mask]
    x_poly  = poly.transform(x_base.reshape(1, -1)).flatten()
    y_pred  = float(np.dot(x_poly, beta))
    se_pred = float(np.sqrt(x_poly @ cov_beta @ x_poly))
    return y_pred, se_pred


def run_predictions(result: dict, base_name: str, genotype_length: int,
                    encoding: str, dropped_features: list,
                    predictions_dir: str) -> None:
    """
    Generate predictions with uncertainty for all 2^L genotypes and save to CSV.
    Output filename mirrors the coefficient base name with '_predictions.csv' suffix.
    """
    os.makedirs(predictions_dir, exist_ok=True)
    output_path = os.path.join(predictions_dir, f"{base_name}_predictions.csv")

    beta     = np.asarray(result["model"].params)
    cov_beta = np.asarray(result["model"].cov_params())
    poly     = result["poly"]

    all_genotypes = generate_all_genotypes(genotype_length)
    total = len(all_genotypes)
    print(f"  Predicting {total:,} genotypes → {output_path}")

    records = []
    for i, geno in enumerate(all_genotypes):
        if (i + 1) % 2000 == 0 or i == total - 1:
            print(f"    {i+1:,}/{total:,} ({100*(i+1)/total:.0f}%)")
        try:
            y_pred, se_pred = predict_with_uncertainty(
                geno, beta, cov_beta, poly, encoding, dropped_features
            )
            records.append({"geno": geno, "predicted_phenotype": y_pred, "prediction_se": se_pred})
        except Exception as e:
            print(f"    WARNING: failed for {geno}: {e}")
            records.append({"geno": geno, "predicted_phenotype": np.nan, "prediction_se": np.nan})

    pred_df = pd.DataFrame(records)
    pred_df.to_csv(output_path, index=False)

    n_ok = pred_df["predicted_phenotype"].notna().sum()
    print(f"  Saved {n_ok:,}/{total:,} predictions  |  "
          f"ŷ mean={pred_df['predicted_phenotype'].mean():.4f}  "
          f"SE mean={pred_df['prediction_se'].mean():.4f}")


def main():
    os.makedirs(MODEL_COEFFS_DIR, exist_ok=True)
    os.makedirs(PREDICTIONS_DIR,  exist_ok=True)

    print(f"Epistasis inference + prediction  |  antigen={ANTIGEN}  encoding={ENCODING}")
    print(f"  Config      : {CONFIG_CSV}")
    print(f"  Coeffs out  : {MODEL_COEFFS_DIR}")
    print(f"  Predict out : {PREDICTIONS_DIR}")
    print("=" * 80)

    # load config settings
    config_df = pd.read_csv(CONFIG_CSV)
    required = {"scope", "phenotype", "order", "predict_phenotypes"}
    missing  = required - set(config_df.columns)
    if missing:
        raise ValueError(f"CONFIG_CSV missing columns: {sorted(missing)}")
    print(f"Config rows: {len(config_df)}")
    print(config_df.to_string(index=False))
    print("=" * 80)

    # Load raw dataset
    raw_df = pd.read_csv(os.path.join(RAW_DATA_DIR, EPISTASIS_CSV))
    raw_df, L = load_and_normalize_genos(raw_df, GENO_COL)
    if POS4_COL not in raw_df.columns:
        raise KeyError(f"'{POS4_COL}' not found in {EPISTASIS_CSV}")
    raw_df[POS4_COL] = raw_df[POS4_COL].astype(int)
    if PHENO_COL not in raw_df.columns:
        raise KeyError(f"'{PHENO_COL}' not found in {EPISTASIS_CSV}")
    print(f"Raw data: rows={len(raw_df)}, L={L}")

    # Load latent datasets
    # Filename: {ANTIGEN}_epistasis_linearized_{link}_v{VERSION}_{ENCODING}.csv
    # VERSION = 1 for scope="full", VERSION = 2 for scope="pos4"
    import glob as _glob
    latent_datasets = {}
    _latent_pattern = os.path.join(LATENT_DATA_DIR, f"{ANTIGEN}_epistasis_linearized_*_v*_{ENCODING}.csv")
    _latent_paths   = sorted(_glob.glob(_latent_pattern))
    print(f"Searching for latent CSVs in: {LATENT_DATA_DIR}")
    print(f"  Pattern : {_latent_pattern}")
    print(f"  Found   : {len(_latent_paths)} file(s)")
    for fpath in _latent_paths:
        base = os.path.basename(fpath)
        stem = (base
                .replace(f"{ANTIGEN}_epistasis_linearized_", "")
                .replace(f"_{ENCODING}.csv", ""))
        df_t, _ = load_and_normalize_genos(pd.read_csv(fpath), GENO_COL)
        if PHENO_COL not in df_t.columns:
            raise KeyError(f"'{PHENO_COL}' not found in {fpath}")
        df_m = raw_df[[GENO_COL, POS4_COL]].merge(
            df_t[[GENO_COL, PHENO_COL]], on=GENO_COL, how="left"
        )
        key = f"latent_{stem}"
        latent_datasets[key] = df_m
        print(f"  Loaded latent: {key}  ({len(df_m)} rows)")

    print(f"Latent datasets loaded: {len(latent_datasets)}")
    print("=" * 80)

    # Process each config row
    for cfg_idx, row in config_df.iterrows():
        scope     = str(row["scope"]).strip()
        phenotype = str(row["phenotype"]).strip()
        max_order = int(row["order"])
        do_predict = str(row["predict_phenotypes"]).strip().upper() == "TRUE"

        print(f"\n[{cfg_idx+1}/{len(config_df)}] scope={scope}  phenotype={phenotype}  "
              f"max_order={max_order}  predict={do_predict}")

        # Select dataset
        if phenotype == "raw":
            df = raw_df.copy()
        else:
            version = 2 if scope == "pos4" else 1
            key     = f"{phenotype}_v{version}"
            if key not in latent_datasets:
                print(f"  WARNING: dataset '{key}' not found — skipping.")
                continue
            df = latent_datasets[key].copy()

        # Scope filter
        if scope == "pos4":
            df = df.loc[df[POS4_COL] == 1].copy()
            print(f"  pos4 filter → {len(df)} rows")
        elif scope == "full":
            print(f"  full dataset → {len(df)} rows")
        else:
            print(f"  WARNING: unknown scope '{scope}' — skipping.")
            continue

        # Drop non-finite phenotypes
        ok     = np.isfinite(df[PHENO_COL].to_numpy(dtype=float))
        n_drop = int(np.sum(~ok))
        if n_drop > 0:
            print(f"  Dropping {n_drop} non-finite phenotype rows")
            df = df.loc[ok].copy()
        if len(df) == 0:
            print("  ERROR: no data after filtering — skipping.")
            continue

        # Encode genotypes
        X_full = encode_genotype_series(df[GENO_COL], ENCODING)
        y      = df[PHENO_COL].to_numpy(dtype=float)

        # For pos4 scope, drop the pos4 feature (constant within subset)
        feature_mask = None
        dropped_features = None
        if scope == "pos4":
            feature_mask = np.ones(X_full.shape[1], dtype=bool)
            feature_mask[POS4_MUT_INDEX] = False
            X = X_full[:, feature_mask]
            feature_labels = [str(i + 1) for i in range(X_full.shape[1]) if feature_mask[i]]
            dropped_features = [POS4_MUT_INDEX]
        else:
            X = X_full
            feature_labels = [str(i + 1) for i in range(X.shape[1])]

        # Fit and export each order 1..max_order
        pheno_clean = phenotype.replace("latent_", "")
        for order in range(1, max_order + 1):
            print(f"\n  Order {order}  (n={len(y)} samples)")

            try:
                result = fit_epistasis_model(X, y, order, feature_labels)
            except Exception as e:
                print(f"  ERROR: fit failed: {e}")
                continue

            base_name = f"{ANTIGEN}_{pheno_clean}_{order}order_{scope}_{ENCODING}"

            # Coefficients
            try:
                export_coefficients(
                    result,
                    os.path.join(MODEL_COEFFS_DIR, f"{base_name}.txt"),
                    alpha=ALPHA,
                )
            except Exception as e:
                print(f"  ERROR: coefficient export failed: {e}")

            # Covariance matrix
            try:
                export_covariance(
                    result,
                    os.path.join(MODEL_COEFFS_DIR, f"{base_name}_covariance.csv"),
                )
            except Exception as e:
                print(f"  ERROR: covariance export failed: {e}")

            # Metadata
            try:
                export_metadata(
                    result,
                    os.path.join(MODEL_COEFFS_DIR, f"{base_name}_metadata.json"),
                    scope, phenotype, order, ENCODING, len(y), feature_mask,
                )
            except Exception as e:
                print(f"  ERROR: metadata export failed: {e}")

            # Predictions (only for the final max_order model, if flagged)
            if do_predict and order == max_order:
                print(f"  Prediction flagged — running for {base_name}")
                try:
                    run_predictions(
                        result, base_name, L, ENCODING, dropped_features,
                        PREDICTIONS_DIR,
                    )
                except Exception as e:
                    print(f"  ERROR: prediction failed: {e}")
            elif do_predict:
                print(f"  Prediction flagged — skipping order {order} (will run at order {max_order})")

    print("\n" + "=" * 80)
    print("Done.")
    print(f"  Model coefficients : {MODEL_COEFFS_DIR}/")
    print(f"  Predictions        : {PREDICTIONS_DIR}/")
    print("=" * 80)


if __name__ == "__main__":
    main()
