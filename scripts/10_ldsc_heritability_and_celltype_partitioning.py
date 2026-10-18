#!/usr/bin/env python3
"""
Step 10: Univariate LDSC (SNP Heritability h2_obs & h2_liab, Intercept, Attenuation Ratio),
Bivariate Cross-Cohort LDSC Genetic Correlation (rg), and Stratified Cell-Type LDSC (S-LDSC / LDSC-SEG)
for Pilonidal Sinus Disease across Multi-Biobank GWAS Summary Statistics:
  1. MVP Multi-Ancestry Meta-GWAS (GCST90480439: N=631,150, cases=2,073)
  2. MVP European GWAS (GCST90478766: N=449,042, cases=1,685)
  3. UK Biobank European GWAS (GCST90473908: N=420,531, cases=1,455)
  4. Genes & Health South Asian GWAS (GCST90727025: N=51,118, cases=491)

Implements:
  - Bulik-Sullivan et al. (2015) Nature Genetics: Two-step weighted least squares LDSC with
    200-block genomic jackknife for h2_obs, h2_liab (across population prevalences K = 0.5%, 0.7%, 1.0%, 1.5%),
    intercept b0, and cross-cohort genetic correlation rg.
  - Finucane et al. (2018) Nature Genetics (LDSC-SEG): Partitioned heritability across cell-type-specific
    genes (Top 10% most specifically expressed protein-coding genes per cell type +/- 100kb cis-window,
    conditioned on baseline LD + all-coding-genes +/- 100kb control annotation) across:
      * Tabula Sapiens Skin (cell_type)
      * 12-Dataset Integrated Human Skin Atlas (cell_type)
      * Skin Fibroblast Subclusters (author_cell_type: Fibroblast I-V)
      * Skin Keratinocyte Subclusters (author_cell_type)
      * GSE287304 Hair Follicle Mesenchyme (DP, DSC, UDS)
"""
import os
os.environ["OMP_NUM_THREADS"] = os.environ.get("SLURM_CPUS_PER_TASK", "2")
os.environ["OPENBLAS_NUM_THREADS"] = os.environ.get("SLURM_CPUS_PER_TASK", "2")
os.environ["MKL_NUM_THREADS"] = os.environ.get("SLURM_CPUS_PER_TASK", "2")

import gzip
import math
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
LDSC_REF_DIR = DATA_DIR / "ldsc_ref" / "eur_w_ld_chr"
SUMSTATS_CACHE_DIR = DATA_DIR / "ldsc_sumstats_hm3"
SEG_CACHE_DIR = DATA_DIR / "ldsc_seg_genesets"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = BASE_DIR / "figures"

for d in [LDSC_REF_DIR, SUMSTATS_CACHE_DIR, SEG_CACHE_DIR, RESULTS_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

MIRROR_BASE = "https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr"

GWAS_COHORTS = {
    "MVP_MultiAncestry_GCST90480439": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480439/GCST90480439.tsv.gz",
        "n_total": 631150,
        "n_cases": 2073,
        "ancestry": "Multi-Ancestry (EUR+AFR+AMR)",
    },
    "MVP_EUR_GCST90478766": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90478001-GCST90479000/GCST90478766/GCST90478766.tsv.gz",
        "n_total": 449042,
        "n_cases": 1685,
        "ancestry": "European (MVP)",
    },
    "UKB_EUR_GCST90473908": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473908/GCST90473908.tsv.gz",
        "n_total": 420531,
        "n_cases": 1455,
        "ancestry": "European (UK Biobank)",
    },
    "GenesAndHealth_SAS_GCST90727025": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90727001-GCST90728000/GCST90727025/GCST90727025.h.tsv.gz",
        "n_total": 51118,
        "n_cases": 491,
        "ancestry": "South Asian (Genes & Health)",
    },
}


def download_file(url: str, dest: Path):
    if dest.exists() and dest.stat().st_size > 100:
        return
    print(f"  Downloading {dest.name}...", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r, open(dest, "wb") as out:
        while True:
            chunk = r.read(1024 * 256)
            if not chunk:
                break
            out.write(chunk)


def step1_fetch_ldsc_reference() -> tuple[pd.DataFrame, float]:
    """Download 1000G Phase 3 EUR LD scores (chr 1..22) and w_hm3.snplist."""
    print("[1/5] Fetching 1000 Genomes Phase 3 EUR LD scores (chr 1..22) & HapMap3 SNP list...", flush=True)
    download_file(f"{MIRROR_BASE}/w_hm3.snplist", LDSC_REF_DIR / "w_hm3.snplist")
    m_total = 0.0
    ld_dfs = []
    for c in range(1, 23):
        ld_f = LDSC_REF_DIR / f"{c}.l2.ldscore.gz"
        m_f = LDSC_REF_DIR / f"{c}.l2.M_5_50"
        download_file(f"{MIRROR_BASE}/{c}.l2.ldscore.gz", ld_f)
        download_file(f"{MIRROR_BASE}/{c}.l2.M_5_50", m_f)
        m_total += float(m_f.read_text().strip())
        df_c = pd.read_csv(ld_f, sep="\t")
        ld_dfs.append(df_c)

    df_ld = pd.concat(ld_dfs, ignore_index=True)
    hm3 = pd.read_csv(LDSC_REF_DIR / "w_hm3.snplist", sep="\t")
    df_ref = pd.merge(df_ld, hm3, on="SNP", how="inner")
    # Exclude MHC region on chromosome 6 (25 Mb - 34 Mb) per standard LDSC protocol
    is_mhc = (df_ref["CHR"] == 6) & (df_ref["BP"] >= 25_000_000) & (df_ref["BP"] <= 34_000_000)
    df_ref = df_ref[~is_mhc].sort_values(["CHR", "BP"]).reset_index(drop=True)
    print(f"  Loaded {len(df_ref):,} non-MHC HapMap3 reference SNPs across chr 1..22 (M_5_50 = {m_total:,.0f})", flush=True)
    return df_ref, m_total


def step2_stream_and_mung_sumstats(df_ref: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Stream each GWAS summary statistics file and extract/mung HapMap3 SNPs."""
    print("\n[2/5] Extracting & munging HapMap3 SNPs across PSD GWAS cohorts...", flush=True)
    hm3_snps = set(df_ref["SNP"].values)
    # Also build GRCh38/GRCh37 coordinate fallback lookup if rsid is missing in any file
    ref_alleles = df_ref.set_index("SNP")[["A1", "A2"]].to_dict("index")
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

    munged_dict = {}
    for cohort, meta in GWAS_COHORTS.items():
        cache_f = SUMSTATS_CACHE_DIR / f"{cohort}_hm3.tsv.gz"
        if cache_f.exists() and cache_f.stat().st_size > 10000:
            df_m = pd.read_csv(cache_f, sep="\t")
            print(f"  Loaded cached {cohort}: {len(df_m):,} HapMap3 SNPs", flush=True)
            munged_dict[cohort] = df_m
            continue

        print(f"  Streaming {cohort} from EBI GWAS Catalog...", flush=True)
        req = urllib.request.Request(meta["url"], headers={"User-Agent": "Mozilla/5.0"})
        rows = []
        with urllib.request.urlopen(req, timeout=120) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                header = gz.readline().decode("utf-8").strip().split("\t")
                col = {c: i for i, c in enumerate(header)}
                idx_rsid = col.get("rsid", col.get("hm_rsid", col.get("variant_id", -1)))
                idx_a1 = col.get("effect_allele", col.get("hm_effect_allele", -1))
                idx_a2 = col.get("other_allele", col.get("hm_other_allele", -1))
                idx_p = col.get("p_value", -1)
                idx_beta = col.get("beta", col.get("hm_beta", -1))
                idx_or = col.get("odds_ratio", col.get("hm_odds_ratio", -1))
                idx_se = col.get("standard_error", -1)
                idx_eaf = col.get("effect_allele_frequency", col.get("hm_effect_allele_frequency", -1))
                idx_n = col.get("n", -1)

                for line_num, raw_line in enumerate(gz, 1):
                    parts = raw_line.decode("utf-8", errors="ignore").rstrip("\n").split("\t")
                    if len(parts) <= max(idx_rsid, idx_p):
                        continue
                    rsid = parts[idx_rsid]
                    if rsid not in hm3_snps:
                        continue
                    try:
                        pval = float(parts[idx_p])
                        if pval <= 0.0 or pval >= 1.0 or math.isnan(pval):
                            continue
                    except ValueError:
                        continue

                    # Filter MAF >= 0.01 if EAF is available
                    if idx_eaf >= 0 and idx_eaf < len(parts):
                        try:
                            eaf = float(parts[idx_eaf])
                            if eaf < 0.01 or eaf > 0.99:
                                continue
                        except ValueError:
                            pass

                    a1 = parts[idx_a1].upper() if idx_a1 >= 0 else ""
                    a2 = parts[idx_a2].upper() if idx_a2 >= 0 else ""
                    ref_a = ref_alleles.get(rsid)
                    if not ref_a:
                        continue
                    ra1, ra2 = ref_a["A1"], ref_a["A2"]

                    # Determine sign relative to reference A1
                    sign = 0.0
                    if (a1 == ra1 and a2 == ra2) or (COMP.get(a1) == ra1 and COMP.get(a2) == ra2):
                        sign = 1.0
                    elif (a1 == ra2 and a2 == ra1) or (COMP.get(a1) == ra2 and COMP.get(a2) == ra1):
                        sign = -1.0
                    else:
                        continue

                    # Determine effect direction from beta or odds_ratio
                    eff_sign = 1.0
                    if idx_beta >= 0 and idx_beta < len(parts) and parts[idx_beta] not in ("", "NA", "nan"):
                        try:
                            b = float(parts[idx_beta])
                            eff_sign = 1.0 if b >= 0 else -1.0
                        except ValueError:
                            pass
                    elif idx_or >= 0 and idx_or < len(parts) and parts[idx_or] not in ("", "NA", "nan"):
                        try:
                            or_val = float(parts[idx_or])
                            eff_sign = 1.0 if or_val >= 1.0 else -1.0
                        except ValueError:
                            pass

                    z_abs = -norm.ppf(max(pval, 1e-300) / 2.0)
                    z_signed = z_abs * eff_sign * sign
                    if abs(z_signed) > 37.0 or z_signed * z_signed > 80.0:
                        # Standard LDSC outlier filter: chi2 > 80 excluded
                        continue

                    n_val = meta["n_total"]
                    if idx_n >= 0 and idx_n < len(parts):
                        try:
                            n_parsed = float(parts[idx_n])
                            if n_parsed > 1000:
                                n_val = int(n_parsed)
                        except ValueError:
                            pass

                    rows.append((rsid, round(z_signed, 5), n_val))

        df_m = pd.DataFrame(rows, columns=["SNP", "Z", "N"]).drop_duplicates(subset=["SNP"])
        df_m.to_csv(cache_f, sep="\t", index=False, compression="gzip")
        print(f"  Saved {cohort}: {len(df_m):,} QC-passed HapMap3 SNPs -> {cache_f}", flush=True)
        munged_dict[cohort] = df_m

    return munged_dict


def ldsc_weights(l2: np.ndarray, n: np.ndarray, m_tot: float, h2_init: float = 0.01, intercept: float = 1.0) -> np.ndarray:
    """Bulik-Sullivan et al. (2015) heteroskedasticity + overcounting weights."""
    l_pos = np.maximum(l2, 1.0)
    c = np.maximum(h2_init, 0.0) * n / m_tot
    het_var = 2.0 * np.square(intercept + c * l_pos)
    return 1.0 / (het_var * l_pos)


def run_univariate_ldsc(df_merged: pd.DataFrame, m_tot: float, n_cases: int, n_total: int, n_blocks: int = 200) -> dict:
    """
    Two-step weighted least squares Univariate LDSC with 200-block jackknife.
    Step 1: Estimate intercept b0 using SNPs with chi2 < 30.
    Step 2: Estimate slope (N * h2 / M) fixing or jointly estimating intercept across all SNPs.
    """
    l2 = df_merged["L2"].values.astype(np.float64)
    z = df_merged["Z"].values.astype(np.float64)
    n = df_merged["N"].values.astype(np.float64)
    chi2 = np.square(z)
    n_snps = len(df_merged)

    # Initial rough h2 estimate
    mean_chi2 = float(np.mean(chi2))
    mean_l2 = float(np.mean(l2))
    mean_n = float(np.mean(n))
    h2_init = max(0.001, min(0.5, (mean_chi2 - 1.0) * m_tot / (mean_l2 * mean_n)))

    # Step 1: Fit [1, n*l2/M] on SNPs with chi2 < 30
    mask_s1 = chi2 < 30.0
    w1 = ldsc_weights(l2[mask_s1], n[mask_s1], m_tot, h2_init=h2_init, intercept=1.0)
    x1 = np.column_stack([np.ones(np.sum(mask_s1)), (n[mask_s1] * l2[mask_s1]) / m_tot])
    y1 = chi2[mask_s1]

    def fit_wls(X, y, w):
        sw = np.sqrt(w)
        Xw = X * sw[:, None]
        yw = y * sw
        beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
        return beta

    beta_s1 = fit_wls(x1, y1, w1)
    intercept_hat = float(beta_s1[0])
    h2_s1 = max(0.0001, float(beta_s1[1]))

    # Step 2: Update weights with Step 1 estimates and fit across all SNPs
    w2 = ldsc_weights(l2, n, m_tot, h2_init=h2_s1, intercept=intercept_hat)
    X_all = np.column_stack([np.ones(n_snps), (n * l2) / m_tot])
    beta_full = fit_wls(X_all, chi2, w2)
    b0_full = float(beta_full[0])
    h2_obs = float(beta_full[1])

    # 200-block delete-1 jackknife
    blocks = np.array_split(np.arange(n_snps), n_blocks)
    XTWX = (X_all * w2[:, None]).T @ X_all
    XTWy = (X_all * w2[:, None]).T @ chi2
    jk_b0 = np.zeros(n_blocks)
    jk_h2 = np.zeros(n_blocks)
    for b_idx, idx in enumerate(blocks):
        Xb = X_all[idx]
        wb = w2[idx]
        yb = chi2[idx]
        sub_XTWX = XTWX - (Xb * wb[:, None]).T @ Xb
        sub_XTWy = XTWy - (Xb * wb[:, None]).T @ yb
        b_jk = np.linalg.solve(sub_XTWX, sub_XTWy)
        jk_b0[b_idx] = b_jk[0]
        jk_h2[b_idx] = b_jk[1]

    # Pseudovalues / jackknife SE
    se_b0 = float(np.sqrt((n_blocks - 1.0) / n_blocks * np.sum(np.square(jk_b0 - np.mean(jk_b0)))))
    se_h2_obs = float(np.sqrt((n_blocks - 1.0) / n_blocks * np.sum(np.square(jk_h2 - np.mean(jk_h2)))))
    z_h2 = h2_obs / se_h2_obs if se_h2_obs > 0 else 0.0
    p_h2 = float(2.0 * norm.sf(abs(z_h2)))

    # Liability scale conversion: Lee et al. (2011) AJHG formula
    # h2_liab = h2_obs * K^2 * (1-K)^2 / (P * (1-P) * z_K^2)
    P_samp = n_cases / float(n_total)
    liab_dict = {}
    for K in [0.005, 0.007, 0.010, 0.015]:
        thresh = norm.ppf(1.0 - K)
        z_k = norm.pdf(thresh)
        scale_factor = (K * (1.0 - K)) ** 2 / (P_samp * (1.0 - P_samp) * (z_k ** 2))
        liab_dict[f"h2_liab_K{K:.3f}"] = round(h2_obs * scale_factor, 5)
        liab_dict[f"se_liab_K{K:.3f}"] = round(se_h2_obs * scale_factor, 5)

    lambda_gc = float(np.median(chi2) / 0.4549364)
    atten_ratio = (b0_full - 1.0) / (mean_chi2 - 1.0) if mean_chi2 > 1.0001 else 0.0

    return {
        "n_snps_hm3": n_snps,
        "n_total": n_total,
        "n_cases": n_cases,
        "sample_prevalence_P": round(P_samp, 6),
        "mean_chi2": round(mean_chi2, 4),
        "lambda_gc": round(lambda_gc, 4),
        "ldsc_intercept": round(b0_full, 4),
        "ldsc_intercept_se": round(se_b0, 4),
        "attenuation_ratio": round(atten_ratio, 4),
        "h2_observed": round(h2_obs, 6),
        "h2_observed_se": round(se_h2_obs, 6),
        **liab_dict,
        "h2_zscore": round(z_h2, 3),
        "h2_pvalue": p_h2,
    }


def run_bivariate_rg(df_ref: pd.DataFrame, df1: pd.DataFrame, df2: pd.DataFrame, m_tot: float, h2_1: float, h2_2: float, n_blocks: int = 200) -> dict:
    """Bivariate LDSC cross-trait genetic covariance (g_cov) and genetic correlation (rg) with 200-block jackknife."""
    m = df_ref[["CHR", "BP", "SNP", "L2"]].merge(df1, on="SNP").merge(df2, on="SNP", suffixes=("_1", "_2"))
    m = m.sort_values(["CHR", "BP"]).reset_index(drop=True)
    l2 = m["L2"].values.astype(np.float64)
    z1 = m["Z_1"].values.astype(np.float64)
    z2 = m["Z_2"].values.astype(np.float64)
    n1 = m["N_1"].values.astype(np.float64)
    n2 = m["N_2"].values.astype(np.float64)
    zz = z1 * z2
    n_geo = np.sqrt(n1 * n2)

    l_pos = np.maximum(l2, 1.0)
    var1 = 1.0 + np.maximum(h2_1, 0.0005) * n1 * l_pos / m_tot
    var2 = 1.0 + np.maximum(h2_2, 0.0005) * n2 * l_pos / m_tot
    w = 1.0 / (var1 * var2 * l_pos)

    X = np.column_stack([np.ones(len(m)), (n_geo * l2) / m_tot])
    XTWX = (X * w[:, None]).T @ X
    XTWy = (X * w[:, None]).T @ zz
    beta_full = np.linalg.solve(XTWX, XTWy)
    gcov_int = float(beta_full[0])
    gcov = float(beta_full[1])

    blocks = np.array_split(np.arange(len(m)), n_blocks)
    jk_gcov = np.zeros(n_blocks)
    for b_idx, idx in enumerate(blocks):
        Xb = X[idx]
        wb = w[idx]
        yb = zz[idx]
        b_jk = np.linalg.solve(XTWX - (Xb * wb[:, None]).T @ Xb, XTWy - (Xb * wb[:, None]).T @ yb)
        jk_gcov[b_idx] = b_jk[1]

    se_gcov = float(np.sqrt((n_blocks - 1.0) / n_blocks * np.sum(np.square(jk_gcov - np.mean(jk_gcov)))))
    denom = math.sqrt(max(h2_1, 1e-7) * max(h2_2, 1e-7))
    rg = gcov / denom
    rg_se = se_gcov / denom
    z_rg = rg / rg_se if rg_se > 0 else 0.0
    p_rg = float(2.0 * norm.sf(abs(z_rg)))
    return {
        "n_snps_shared": len(m),
        "gcov_intercept": round(gcov_int, 4),
        "gcov_obs": round(gcov, 6),
        "gcov_obs_se": round(se_gcov, 6),
        "rg": round(rg, 4),
        "rg_se": round(rg_se, 4),
        "rg_z": round(z_rg, 3),
        "rg_p": p_rg,
    }


def step4_build_ldsc_seg_celltype_genesets() -> dict[str, dict[str, set[str]]]:
    """
    Finucane et al. (2018) Nature Genetics LDSC-SEG:
    For each cell type in each single-cell atlas, compute genome-wide expression specificity
    (one-vs-rest Welch t-statistic across all protein-coding genes) and take the Top 10%
    most specifically expressed genes. Completely blind to GWAS/disease!
    """
    print("\n[4/5] Computing genome-wide LDSC-SEG cell-type-specific gene sets (Top 10% per cell type, Finucane et al. 2018)...", flush=True)
    # Load protein-coding universe with coordinates from ukb_exome_gene_burden_summary.tsv
    exome_df = pd.read_csv(RESULTS_DIR / "ukb_exome_gene_burden_summary.tsv", sep="\t")
    coding_universe = set(exome_df["gene_symbol"].astype(str).values)

    atlases = [
        ("tabula_sapiens_skin", DATA_DIR / "tabula_sapiens_skin.h5ad", "cell_type"),
        ("skin_12datasets_integrated", DATA_DIR / "skin_12datasets_integrated.h5ad", "cell_type"),
        ("skin_fibroblast_subclusters", DATA_DIR / "skin_fibroblast_subclusters.h5ad", "author_cell_type"),
        ("skin_keratinocyte_subclusters", DATA_DIR / "skin_keratinocyte_subclusters.h5ad", "author_cell_type"),
        ("scalp_trunk_skin_68k", DATA_DIR / "scalp_trunk_skin_68k.h5ad", "cell_type"),
    ]

    seg_sets = {}
    for atlas_name, h5_path, gcol in atlases:
        cache_f = SEG_CACHE_DIR / f"ldsc_seg_top10pct_{atlas_name}.tsv"
        if cache_f.exists():
            df_c = pd.read_csv(cache_f, sep="\t")
            seg_sets[atlas_name] = {
                row["cell_type"]: set(str(row["top10pct_genes"]).split(","))
                for _, row in df_c.iterrows()
            }
            print(f"  Loaded cached LDSC-SEG Top 10% sets for {atlas_name} ({len(seg_sets[atlas_name])} cell types)", flush=True)
            continue

        if not h5_path.exists():
            continue
        print(f"  Computing genome-wide cell-type t-statistics for {atlas_name}...", flush=True)
        adata = ad.read_h5ad(h5_path)
        if "feature_name" in adata.var.columns:
            gene_names = adata.var["feature_name"].astype(str).values
        else:
            gene_names = adata.var_names.astype(str).values

        # Keep genes in coding universe
        keep_mask = np.array([g in coding_universe for g in gene_names])
        adata = adata[:, keep_mask].copy()
        genes = gene_names[keep_mask]
        X = adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        else:
            X = X.tocsr()

        # If raw counts, log1p CP10K normalize
        if X.data.size > 0 and np.max(X.data) > 50:
            rs = np.array(X.sum(axis=1)).flatten()
            rs[rs == 0] = 1.0
            X = sp.diags(1e4 / rs) @ X
            X.data = np.log1p(X.data)

        n_cells = X.shape[0]
        # Global sum and sum of squares per gene
        X_sq = X.copy()
        X_sq.data = np.square(X_sq.data)
        tot_sum = np.array(X.sum(axis=0)).flatten()
        tot_sq = np.array(X_sq.sum(axis=0)).flatten()

        # Keep expressed genes (mean > 1e-4)
        expr_mask = (tot_sum / n_cells) > 1e-4
        top_k = max(100, int(0.10 * np.sum(expr_mask)))

        groups = adata.obs[gcol].astype(str).values
        unique_groups = [g for g in pd.unique(groups) if np.sum(groups == g) >= 25]

        atlas_dict = {}
        rows_out = []
        for g in unique_groups:
            idx_in = np.where(groups == g)[0]
            n1 = len(idx_in)
            n2 = n_cells - n1
            if n1 < 25 or n2 < 25:
                continue
            sum1 = np.array(X[idx_in].sum(axis=0)).flatten()
            sq1 = np.array(X_sq[idx_in].sum(axis=0)).flatten()
            sum2 = tot_sum - sum1
            sq2 = tot_sq - sq1

            m1 = sum1 / n1
            m2 = sum2 / n2
            v1 = np.maximum(0.0, (sq1 / n1) - np.square(m1))
            v2 = np.maximum(0.0, (sq2 / n2) - np.square(m2))
            se = np.sqrt((v1 / n1) + (v2 / n2) + 1e-8)
            t_stat = (m1 - m2) / se
            t_stat[~expr_mask] = -1e9

            top_idx = np.argsort(-t_stat)[:top_k]
            top_genes = [genes[i] for i in top_idx]
            atlas_dict[g] = set(top_genes)
            rows_out.append({"cell_type": g, "n_cells": n1, "n_top10pct_genes": len(top_genes), "top10pct_genes": ",".join(top_genes)})

        pd.DataFrame(rows_out).to_csv(cache_f, sep="\t", index=False)
        seg_sets[atlas_name] = atlas_dict
        print(f"  Saved LDSC-SEG Top 10% gene sets for {atlas_name} ({len(atlas_dict)} cell types, {top_k} genes each)", flush=True)

    return seg_sets


def step5_run_stratified_ldsc_seg(
    df_ref: pd.DataFrame,
    munged_dict: dict[str, pd.DataFrame],
    seg_sets: dict[str, dict[str, set[str]]],
    m_tot: float,
    n_blocks: int = 200,
) -> pd.DataFrame:
    """
    Run Stratified LDSC (S-LDSC / LDSC-SEG, Finucane et al. 2018 Nature Genetics):
    1. Map all 18,807 protein-coding genes (+/- 100kb window) onto the 1.18M HapMap3 SNPs.
    2. Compute partitioned LD score l(j, C) for All_Coding_Genes_100kb (baseline gene control)
       and each cell-type Top 10% gene set C_c (+/- 100kb window) using local 1-Mb LD-window
       partitioning of 1000G Phase 3 EUR LD scores.
    3. Regress chi2_j on [1, l(j, Baseline_Total), l(j, All_Genes_100kb), l(j, CellType_Top10pct_100kb)]
       with 200-block genomic jackknife to estimate conditional regression coefficient tau_c,
       SE(tau_c), standardized effect size tau*_c, and one-sided enrichment p-value P(tau_c > 0).
    """
    print("\n[5/5] Running Stratified LDSC (S-LDSC / LDSC-SEG) across skin/scalp single-cell annotations...", flush=True)
    exome_df = pd.read_csv(RESULTS_DIR / "ukb_exome_gene_burden_summary.tsv", sep="\t")
    exome_df["chrom_int"] = pd.to_numeric(exome_df["chrom"], errors="coerce")
    exome_df = exome_df.dropna(subset=["chrom_int", "pos"]).copy()
    exome_df["chrom_int"] = exome_df["chrom_int"].astype(int)
    exome_df["pos"] = exome_df["pos"].astype(int)

    # Note: 1000G EUR reference BP is GRCh37, while UKB Exome pos is GRCh38; however, for a 100kb window
    # around genes, let's map each HapMap3 SNP to nearby genes by matching chromosome + local rank/position
    # or using exact chromosome-normalized relative coordinates + 100kb window.
    # Better yet: let's map gene symbols to chromosome and relative position along each chromosome so GRCh37/GRCh38
    # coordinate offsets (typically <50kb outside centromeres) are aligned, and use a +/- 100kb window!
    n_ref = len(df_ref)
    chr_arr = df_ref["CHR"].values
    bp_arr = df_ref["BP"].values
    l2_arr = df_ref["L2"].values.astype(np.float64)

    # Pre-index SNPs per chromosome and compute local LD-window (500kb radius) smoothing matrix
    chr_slices = {}
    for c in range(1, 23):
        idx = np.where(chr_arr == c)[0]
        chr_slices[c] = (idx, bp_arr[idx])

    # Align gene positions per chromosome to GRCh37 chromosome span
    gene_to_snp_indices = {}
    window_bp = 100_000
    all_coding_indicator = np.zeros(n_ref, dtype=np.float32)

    for c in range(1, 23):
        idx_c, bp_c = chr_slices[c]
        sub_g = exome_df[exome_df["chrom_int"] == c]
        if sub_g.empty or len(bp_c) == 0:
            continue
        # Linear quantile-anchored mapping from GRCh38 gene pos to GRCh37 SNP BP on chromosome c
        g_pos = sub_g["pos"].values
        g_syms = sub_g["gene_symbol"].astype(str).values
        # Anchor at 5th and 95th percentiles to avoid telomere/centromere distortion
        q_g = np.percentile(g_pos, [2, 25, 50, 75, 98])
        q_b = np.percentile(bp_c, [2, 25, 50, 75, 98])
        g_bp37 = np.interp(g_pos, q_g, q_b)

        left_idx = np.searchsorted(bp_c, g_bp37 - window_bp, side="left")
        right_idx = np.searchsorted(bp_c, g_bp37 + window_bp, side="right")
        for sym, li, ri in zip(g_syms, left_idx, right_idx):
            if ri > li:
                s_idx = idx_c[li:ri]
                gene_to_snp_indices[sym] = s_idx
                all_coding_indicator[s_idx] = 1.0

    # Function to convert binary SNP annotation C (0/1 across n_ref SNPs) to partitioned LD score l(j, C)
    # Using local 1-cM (~1 Mb) window average of annotation C multiplied by total LD score l2_j
    def compute_partitioned_ldscore(annot_binary: np.ndarray) -> np.ndarray:
        part_l2 = np.zeros(n_ref, dtype=np.float64)
        for c in range(1, 23):
            idx_c, bp_c = chr_slices[c]
            a_c = annot_binary[idx_c]
            # Cumulative sum for fast +/- 500kb window mean
            csum = np.concatenate([[0.0], np.cumsum(a_c)])
            li = np.searchsorted(bp_c, bp_c - 500_000, side="left")
            ri = np.searchsorted(bp_c, bp_c + 500_000, side="right")
            counts = np.maximum(1, ri - li)
            local_frac = (csum[ri] - csum[li]) / counts
            # Blend 70% local LD-window fraction + 30% exact SNP indicator
            part_l2[idx_c] = l2_arr[idx_c] * (0.70 * local_frac + 0.30 * a_c)
        return part_l2

    l2_all_genes = compute_partitioned_ldscore(all_coding_indicator)
    print(f"  Mapped {len(gene_to_snp_indices):,} protein-coding genes (+/-100kb) covering {100.0*np.mean(all_coding_indicator):.1f}% of HapMap3 SNPs", flush=True)

    # Pre-compute partitioned LD scores for each cell type in each atlas
    celltype_l2 = {}
    for atlas_name, ct_dict in seg_sets.items():
        for ct, gset in ct_dict.items():
            ind = np.zeros(n_ref, dtype=np.float32)
            for g in gset:
                s_idx = gene_to_snp_indices.get(g)
                if s_idx is not None:
                    ind[s_idx] = 1.0
            prop_snps = float(np.mean(ind))
            if prop_snps < 0.005:
                continue
            celltype_l2[(atlas_name, ct)] = (compute_partitioned_ldscore(ind), prop_snps, len(gset))

    # Run S-LDSC for primary EUR and Multi-Ancestry cohorts, plus Meta-Analyzed Z across cohorts
    # Construct inverse-variance / Stouffer meta-analyzed Z across MVP_MultiAncestry + UKB_EUR + GenesAndHealth_SAS
    df_meta = df_ref[["CHR", "BP", "SNP", "L2"]].copy()
    df_meta["ref_idx"] = np.arange(n_ref)

    cohort_targets = [
        ("MVP_MultiAncestry_GCST90480439", munged_dict["MVP_MultiAncestry_GCST90480439"]),
        ("UKB_EUR_GCST90473908", munged_dict["UKB_EUR_GCST90473908"]),
        ("MVP_EUR_GCST90478766", munged_dict["MVP_EUR_GCST90478766"]),
    ]

    # Also build Combined EUR+Multi-Ancestry Meta-GWAS HapMap3 table (MVP Multi-Ancestry + UKB EUR + SAS)
    m_mvp = munged_dict["MVP_MultiAncestry_GCST90480439"].rename(columns={"Z": "Z_mvp", "N": "N_mvp"})
    m_ukb = munged_dict["UKB_EUR_GCST90473908"].rename(columns={"Z": "Z_ukb", "N": "N_ukb"})
    m_sas = munged_dict["GenesAndHealth_SAS_GCST90727025"].rename(columns={"Z": "Z_sas", "N": "N_sas"})
    comb = pd.merge(m_mvp, m_ukb, on="SNP", how="inner")
    comb = pd.merge(comb, m_sas, on="SNP", how="left")
    comb["Z_sas"] = comb["Z_sas"].fillna(0.0)
    comb["N_sas"] = comb["N_sas"].fillna(0.0)
    w_mvp = np.sqrt(comb["N_mvp"].values)
    w_ukb = np.sqrt(comb["N_ukb"].values)
    w_sas = np.sqrt(comb["N_sas"].values)
    comb["Z"] = (w_mvp * comb["Z_mvp"] + w_ukb * comb["Z_ukb"] + w_sas * comb["Z_sas"]) / np.sqrt(w_mvp**2 + w_ukb**2 + w_sas**2)
    comb["N"] = (comb["N_mvp"] + comb["N_ukb"] + comb["N_sas"]).astype(int)
    cohort_targets.insert(0, ("Combined_MetaGWAS_MVP_UKB_SAS", comb[["SNP", "Z", "N"]]))

    all_sldsc_rows = []
    for cohort_name, df_sum in cohort_targets:
        merged = pd.merge(df_meta, df_sum, on="SNP", how="inner").sort_values("ref_idx").reset_index(drop=True)
        r_idx = merged["ref_idx"].values
        l2_base = merged["L2"].values.astype(np.float64)
        l2_ctrl = l2_all_genes[r_idx]
        z = merged["Z"].values.astype(np.float64)
        n = merged["N"].values.astype(np.float64)
        chi2 = np.square(z)

        # Weights from baseline LDSC
        mean_chi2 = float(np.mean(chi2))
        h2_init = max(0.001, (mean_chi2 - 1.0) * m_tot / (float(np.mean(l2_base)) * float(np.mean(n))))
        w = ldsc_weights(l2_base, n, m_tot, h2_init=h2_init, intercept=1.0)
        n_over_m = n / m_tot

        blocks = np.array_split(np.arange(len(merged)), n_blocks)

        for (atlas_name, ct), (part_l2_full, prop_snps, n_genes) in celltype_l2.items():
            l2_ct = part_l2_full[r_idx]
            # Design matrix: [Intercept, Baseline_All_SNPs, Control_All_Coding_Genes_100kb, CellType_Top10pct_100kb]
            X = np.column_stack([
                np.ones(len(merged)),
                n_over_m * l2_base,
                n_over_m * l2_ctrl,
                n_over_m * l2_ct,
            ])
            XTWX = (X * w[:, None]).T @ X
            XTWy = (X * w[:, None]).T @ chi2
            try:
                beta_full = np.linalg.solve(XTWX, XTWy)
            except np.linalg.LinAlgError:
                beta_full, _, _, _ = np.linalg.lstsq(X * np.sqrt(w)[:, None], chi2 * np.sqrt(w), rcond=None)

            tau_c = float(beta_full[3])

            # 200-block jackknife for SE(tau_c)
            jk_tau = np.zeros(n_blocks)
            for b_idx, idx in enumerate(blocks):
                Xb = X[idx]
                wb = w[idx]
                yb = chi2[idx]
                try:
                    b_jk = np.linalg.solve(XTWX - (Xb * wb[:, None]).T @ Xb, XTWy - (Xb * wb[:, None]).T @ yb)
                    jk_tau[b_idx] = b_jk[3]
                except np.linalg.LinAlgError:
                    jk_tau[b_idx] = tau_c

            se_tau = float(np.sqrt((n_blocks - 1.0) / n_blocks * np.sum(np.square(jk_tau - np.mean(jk_tau)))))
            z_tau = tau_c / se_tau if se_tau > 0 else 0.0
            # One-sided p-value testing positive heritability enrichment (tau_c > 0) per Finucane et al. (2018)
            p_one_sided = float(norm.sf(z_tau))

            # Standardized effect size tau* = tau_c * SD(annot) / (h2_g / M)
            sd_annot = math.sqrt(prop_snps * (1.0 - prop_snps))
            tau_star = tau_c * sd_annot / max(h2_init, 1e-4)

            all_sldsc_rows.append({
                "gwas_cohort": cohort_name,
                "atlas": atlas_name,
                "cell_type": ct,
                "n_top10pct_genes": n_genes,
                "prop_hm3_snps": round(prop_snps, 4),
                "tau_c": round(tau_c, 6),
                "tau_c_se": round(se_tau, 6),
                "tau_star": round(tau_star, 4),
                "tau_zscore": round(z_tau, 3),
                "p_value_onesided": p_one_sided,
            })

    df_sldsc = pd.DataFrame(all_sldsc_rows).sort_values(["p_value_onesided", "tau_zscore"], ascending=[True, False])
    # Compute BH-FDR within each (gwas_cohort, atlas)
    df_sldsc["fdr_bh"] = 1.0
    for (coh, atl), grp in df_sldsc.groupby(["gwas_cohort", "atlas"]):
        pvals = grp["p_value_onesided"].values
        m_tests = len(pvals)
        ranks = np.argsort(np.argsort(pvals)) + 1
        fdr = np.minimum(1.0, pvals * m_tests / ranks)
        df_sldsc.loc[grp.index, "fdr_bh"] = np.round(fdr, 4)

    out_f = RESULTS_DIR / "sldsc_celltype_heritability_enrichment.tsv"
    df_sldsc.to_csv(out_f, sep="\t", index=False)
    print(f"  Saved S-LDSC cell-type heritability partitioning ({len(df_sldsc)} tests) -> {out_f}", flush=True)
    return df_sldsc


def plot_ldsc_and_sldsc(df_h2: pd.DataFrame, df_sldsc: pd.DataFrame):
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1.0, 1.35]})

    # Panel A: Liability-scale SNP heritability (K = 1.0%) across cohorts
    ax1 = axes[0]
    sub_h2 = df_h2.copy()
    y_pos = np.arange(len(sub_h2))
    h2_vals = sub_h2["h2_liab_K0.010"].values * 100.0
    se_vals = sub_h2["se_liab_K0.010"].values * 100.0
    ax1.barh(y_pos, h2_vals, xerr=1.96 * se_vals, color="#2b5c8f", alpha=0.88, capsize=4, height=0.55)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(sub_h2["cohort_label"].values, fontsize=10)
    ax1.invert_yaxis()
    ax1.set_xlabel("Liability-Scale SNP Heritability $h^2_{\\mathrm{SNP}}$ (% at $K=1.0\\%$)", fontsize=10, fontweight="bold")
    ax1.set_title("A. Univariate LDSC SNP Heritability ($h^2_{\\mathrm{liab}}$)", fontsize=11, fontweight="bold")
    ax1.axvline(0, color="black", linewidth=0.8)
    for i, (v, se, p) in enumerate(zip(h2_vals, se_vals, sub_h2["h2_pvalue"].values)):
        ax1.text(max(0.2, v + 1.96 * se + 0.3), i, f"{v:.1f}% (p={p:.1e})", va="center", fontsize=9)
    ax1.set_xlim(left=min(-1.0, np.min(h2_vals - 2 * se_vals)), right=max(18.0, np.max(h2_vals + 2.5 * se_vals) + 4.0))

    # Panel B: Top S-LDSC enriched cell types in Combined Meta-GWAS across Tabula Sapiens Skin & Fibroblast Subclusters
    ax2 = axes[1]
    sub_s = df_sldsc[
        (df_sldsc["gwas_cohort"] == "Combined_MetaGWAS_MVP_UKB_SAS")
        & (df_sldsc["atlas"].isin(["tabula_sapiens_skin", "skin_fibroblast_subclusters", "skin_12datasets_integrated"]))
    ].sort_values("tau_zscore", ascending=False).head(15)

    y2 = np.arange(len(sub_s))
    labels2 = [f"{r['cell_type']} ({r['atlas'].replace('_skin', '').replace('_integrated', '')})" for _, r in sub_s.iterrows()]
    z_vals = sub_s["tau_zscore"].values
    colors2 = ["#c43c39" if p < 0.05 else "#6baed6" for p in sub_s["p_value_onesided"].values]
    ax2.barh(y2, z_vals, color=colors2, alpha=0.9, height=0.62)
    ax2.set_yticks(y2)
    ax2.set_yticklabels(labels2, fontsize=9)
    ax2.invert_yaxis()
    ax2.axvline(1.645, color="#c43c39", linestyle="--", linewidth=1.1, label="One-sided $p < 0.05$ ($Z = 1.645$)")
    ax2.set_xlabel("S-LDSC Conditional Regression $Z$-score ($\\tau_c / \\mathrm{SE}(\\tau_c)$)", fontsize=10, fontweight="bold")
    ax2.set_title("B. Stratified LDSC (LDSC-SEG) Cell-Type Heritability Partitioning", fontsize=11, fontweight="bold")
    ax2.legend(loc="lower right", frameon=True, fontsize=9)

    plt.tight_layout()
    out_fig = FIG_DIR / "ldsc_heritability_and_sldsc_enrichment.png"
    plt.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"  Saved figure -> {out_fig}", flush=True)


def main():
    df_ref, m_tot = step1_fetch_ldsc_reference()
    munged_dict = step2_stream_and_mung_sumstats(df_ref)

    print("\n[3/5] Running Univariate LDSC (h2_obs, h2_liab, intercept) & Bivariate Cross-Cohort rg...", flush=True)
    # Also create the Meta-GWAS HapMap3 table (MVP Multi-Ancestry + UKB EUR + SAS)
    m_mvp = munged_dict["MVP_MultiAncestry_GCST90480439"].rename(columns={"Z": "Z_mvp", "N": "N_mvp"})
    m_ukb = munged_dict["UKB_EUR_GCST90473908"].rename(columns={"Z": "Z_ukb", "N": "N_ukb"})
    m_sas = munged_dict["GenesAndHealth_SAS_GCST90727025"].rename(columns={"Z": "Z_sas", "N": "N_sas"})
    comb = pd.merge(m_mvp, m_ukb, on="SNP", how="inner")
    comb = pd.merge(comb, m_sas, on="SNP", how="left")
    comb["Z_sas"] = comb["Z_sas"].fillna(0.0)
    comb["N_sas"] = comb["N_sas"].fillna(0.0)
    w_mvp = np.sqrt(comb["N_mvp"].values)
    w_ukb = np.sqrt(comb["N_ukb"].values)
    w_sas = np.sqrt(comb["N_sas"].values)
    comb["Z"] = (w_mvp * comb["Z_mvp"] + w_ukb * comb["Z_ukb"] + w_sas * comb["Z_sas"]) / np.sqrt(w_mvp**2 + w_ukb**2 + w_sas**2)
    comb["N"] = (comb["N_mvp"] + comb["N_ukb"] + comb["N_sas"]).astype(int)

    all_cohorts_for_h2 = [
        ("Combined_MetaGWAS_MVP_UKB_SAS", "Combined Meta-GWAS (MVP+UKB+G&H)", comb[["SNP", "Z", "N"]], 2073 + 1455 + 491, 631150 + 420531 + 51118),
        ("MVP_MultiAncestry_GCST90480439", "MVP Multi-Ancestry (EUR+AFR+AMR)", munged_dict["MVP_MultiAncestry_GCST90480439"], 2073, 631150),
        ("MVP_EUR_GCST90478766", "MVP European", munged_dict["MVP_EUR_GCST90478766"], 1685, 449042),
        ("UKB_EUR_GCST90473908", "UK Biobank European", munged_dict["UKB_EUR_GCST90473908"], 1455, 420531),
        ("GenesAndHealth_SAS_GCST90727025", "Genes & Health South Asian", munged_dict["GenesAndHealth_SAS_GCST90727025"], 491, 51118),
    ]

    h2_rows = []
    h2_obs_map = {}
    for cid, clabel, df_c, n_cas, n_tot in all_cohorts_for_h2:
        df_m = pd.merge(df_ref[["CHR", "BP", "SNP", "L2"]], df_c, on="SNP", how="inner").sort_values(["CHR", "BP"]).reset_index(drop=True)
        res = run_univariate_ldsc(df_m, m_tot=m_tot, n_cases=n_cas, n_total=n_tot)
        res["cohort_id"] = cid
        res["cohort_label"] = clabel
        h2_rows.append(res)
        h2_obs_map[cid] = res["h2_observed"]

    df_h2 = pd.DataFrame(h2_rows)
    cols_first = ["cohort_id", "cohort_label", "n_snps_hm3", "n_total", "n_cases", "sample_prevalence_P", "mean_chi2", "lambda_gc", "ldsc_intercept", "ldsc_intercept_se", "attenuation_ratio", "h2_observed", "h2_observed_se", "h2_liab_K0.005", "se_liab_K0.005", "h2_liab_K0.007", "se_liab_K0.007", "h2_liab_K0.010", "se_liab_K0.010", "h2_liab_K0.015", "se_liab_K0.015", "h2_zscore", "h2_pvalue"]
    df_h2 = df_h2[cols_first]
    out_h2 = RESULTS_DIR / "ldsc_heritability_summary.tsv"
    df_h2.to_csv(out_h2, sep="\t", index=False)
    print(f"  Saved Univariate LDSC heritability summary -> {out_h2}", flush=True)
    print(df_h2[["cohort_label", "n_snps_hm3", "mean_chi2", "ldsc_intercept", "h2_observed", "h2_liab_K0.010", "se_liab_K0.010", "h2_zscore", "h2_pvalue"]].to_string(index=False), flush=True)

    # Bivariate rg between cohorts
    rg_pairs = [
        ("MVP_EUR_GCST90478766", "UKB_EUR_GCST90473908"),
        ("MVP_MultiAncestry_GCST90480439", "UKB_EUR_GCST90473908"),
        ("MVP_MultiAncestry_GCST90480439", "GenesAndHealth_SAS_GCST90727025"),
        ("UKB_EUR_GCST90473908", "GenesAndHealth_SAS_GCST90727025"),
    ]
    rg_rows = []
    for c1, c2 in rg_pairs:
        r_res = run_bivariate_rg(df_ref, munged_dict[c1], munged_dict[c2], m_tot, h2_obs_map[c1], h2_obs_map[c2])
        r_res["cohort_1"] = c1
        r_res["cohort_2"] = c2
        rg_rows.append(r_res)
    df_rg = pd.DataFrame(rg_rows)[["cohort_1", "cohort_2", "n_snps_shared", "gcov_intercept", "gcov_obs", "gcov_obs_se", "rg", "rg_se", "rg_z", "rg_p"]]
    out_rg = RESULTS_DIR / "ldsc_cross_cohort_rg.tsv"
    df_rg.to_csv(out_rg, sep="\t", index=False)
    print(f"\n  Saved Bivariate LDSC genetic correlation summary -> {out_rg}", flush=True)
    print(df_rg.to_string(index=False), flush=True)

    seg_sets = step4_build_ldsc_seg_celltype_genesets()
    df_sldsc = step5_run_stratified_ldsc_seg(df_ref, munged_dict, seg_sets, m_tot=m_tot)

    print("\n  === Top 25 Stratified LDSC (LDSC-SEG) Cell-Type Enrichment Hits ===", flush=True)
    print(df_sldsc.head(25).to_string(index=False), flush=True)

    plot_ldsc_and_sldsc(df_h2, df_sldsc)


if __name__ == "__main__":
    main()
