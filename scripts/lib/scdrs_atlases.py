#!/usr/bin/env python3
"""
Helper module (formerly step 3). Downloads human skin & hair-follicle single-cell atlases (CELLxGENE Tabula Sapiens Skin,
Scalp/Trunk Skin Atlas, Keratinocyte & Fibroblast subcluster atlases, and GEO Hair Follicle
Mesenchymal scRNA-seq GSE287304), run scDRS (1,000 Monte Carlo control sets + cell-type association
and heterogeneity tests), and compute gene-by-cell-type specificity matrices.

Current role: only its atlas helpers (ATLAS_URLS, prepare_adata_with_symbols, choose_group_columns, load_gs_file) are used,
imported by step 03b. Its own main run (early MVP-only gene sets) is superseded and its outputs were removed.
"""
import os
import urllib.request
import warnings
from pathlib import Path
import numpy as np
if not hasattr(np, "float_"):
    np.float_ = np.float64
import scipy.stats
_orig_rankdata = scipy.stats.rankdata
def _int_rankdata(a, method="average", *args, **kwargs):
    res = _orig_rankdata(np.nan_to_num(a, nan=0.0), method=method, *args, **kwargs)
    return res.astype(np.intp) if method == "ordinal" else res
scipy.stats.rankdata = _int_rankdata
import pandas as pd
import scanpy as sc
import anndata as ad
import scdrs

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parents[2]   # scripts/lib/ -> repo root
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = BASE_DIR / "figures"
DATA_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

ATLAS_URLS = {
    "tabula_sapiens_skin": (
        "https://datasets.cellxgene.cziscience.com/c240d5bc-a521-4330-a76b-7d566e939d87.h5ad",
        "Tabula Sapiens - Skin (17,786 cells)"
    ),
    "skin_12datasets_integrated": (
        "https://datasets.cellxgene.cziscience.com/19292dfa-1aa8-43be-966f-4172977141e0.h5ad",
        "Human Skin 12-Dataset Atlas (15,243 cells: Dermis + Epidermis)"
    ),
    "skin_keratinocyte_subclusters": (
        "https://datasets.cellxgene.cziscience.com/9662b6de-67f5-4698-a0e6-e376787a7b68.h5ad",
        "Human Skin Keratinocyte Subclusters (9,471 cells)"
    ),
    "skin_fibroblast_subclusters": (
        "https://datasets.cellxgene.cziscience.com/99c00c97-629f-49a4-9b3b-46a930a6d09d.h5ad",
        "Human Skin Fibroblast Subclusters (2,303 cells)"
    ),
    "scalp_trunk_skin_68k": (
        "https://datasets.cellxgene.cziscience.com/688f626d-bb26-40dd-8465-83022a5b2c5f.h5ad",
        "Human Scalp, Trunk & Foreskin Atlas (68,036 cells)"
    ),
}


def load_gs_file(gs_path: Path):
    gs_dict = {}
    with open(gs_path) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            trait, gstr = line.strip().split("\t")
            pairs = [x.split(":") for x in gstr.split(",") if ":" in x]
            gs_dict[trait] = ([p[0] for p in pairs], [float(p[1]) for p in pairs])
    return gs_dict


def prepare_adata_with_symbols(h5ad_path: Path) -> ad.AnnData:
    adata_orig = sc.read_h5ad(h5ad_path)
    obsp_saved = dict(adata_orig.obsp)
    obsm_saved = dict(adata_orig.obsm)
    # Use raw counts if available in .raw, then normalize & log1p for scDRS
    if adata_orig.raw is not None:
        adata = adata_orig.raw.to_adata()
        for k, v in obsp_saved.items():
            adata.obsp[k] = v
        for k, v in obsm_saved.items():
            adata.obsm[k] = v
    else:
        adata = adata_orig
    # Set var_names to gene symbols (feature_name in CELLxGENE h5ad)
    if "feature_name" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names.astype(str)
        adata.var_names = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()
    # Ensure standard log-normalized counts
    x_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
    if x_max > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if "connectivities" not in adata.obsp:
        use_rep = None
        for rep in ["X_pca", "X_scVI", "X_harmony", "X_umap"]:
            if rep in adata.obsm:
                use_rep = rep
                break
        if use_rep is None:
            sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))
            use_rep = "X_pca"
        sc.pp.neighbors(adata, use_rep=use_rep)
    return adata


def choose_group_columns(adata: ad.AnnData):
    cols = []
    for candidate in ["cell_type", "author_cell_type", "cluster", "subcluster", "cell_type_fine", "tissue"]:
        if candidate in adata.obs.columns and adata.obs[candidate].nunique() > 1:
            cols.append(candidate)
    if not cols:
        # fallback to any categorical column with 2..60 levels
        for c in adata.obs.columns:
            n_u = adata.obs[c].nunique()
            if 2 <= n_u <= 60 and c not in ["donor_id", "sample", "batch", "assay"]:
                cols.append(c)
                break
    return cols


def run_scdrs_on_atlas(atlas_key: str, url: str, desc: str, gs_dict: dict):
    h5ad_path = DATA_DIR / f"{atlas_key}.h5ad"
    if not h5ad_path.exists():
        print(f"\n[Download] Fetching {desc} -> {h5ad_path.name} ...")
        urllib.request.urlretrieve(url, h5ad_path)

    print(f"\n[scDRS] Loading {desc} ({h5ad_path.name})...")
    adata = prepare_adata_with_symbols(h5ad_path)
    print(f"  Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes | obs cols: {list(adata.obs.columns)}")

    group_cols = choose_group_columns(adata)
    print(f"  Grouping columns for downstream association: {group_cols}")

    # Preprocess for scDRS (computes gene-level mean/var bins for matched control gene sets)
    scdrs.preprocess(adata, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)

    all_group_stats = []
    cell_scores_summary = {}

    for trait, (gene_list, gene_weights) in gs_dict.items():
        present = [(g, w) for g, w in zip(gene_list, gene_weights) if g in adata.var_names]
        if len(present) < 3:
            print(f"  Skipping {trait}: only {len(present)} genes present")
            continue
        g_list = [p[0] for p in present]
        w_list = [p[1] for p in present]
        print(f"  Scoring trait '{trait}' ({len(g_list)}/{len(gene_list)} genes matched)...")

        df_score = scdrs.score_cell(
            data=adata,
            gene_list=g_list,
            gene_weight=w_list,
            ctrl_match_key="mean_var",
            n_ctrl=1000,
            weight_opt="vs",
            return_ctrl_raw_score=False,
            return_ctrl_norm_score=True,
            verbose=False,
        )
        df_score = df_score.fillna(0.0)

        adata.obs[f"scdrs_norm_{trait}"] = df_score["norm_score"].values
        adata.obs[f"scdrs_pval_{trait}"] = df_score["mc_pval"].values
        fdr = scdrs.util.multipletests(df_score["mc_pval"].values, method="fdr_bh")[1]
        adata.obs[f"scdrs_fdr_{trait}"] = fdr

        for gcol in group_cols:
            # Perform official scDRS downstream cell-type group analysis (assoc_mcp & hetero_mcp)
            dict_df = scdrs.method.downstream_group_analysis(
                adata=adata,
                df_full_score=df_score,
                group_cols=[gcol],
                fdr_thresholds=[0.05, 0.1, 0.2],
            )
            df_res = dict_df[gcol].reset_index().rename(columns={"group": "group_level", gcol: "group_level", "index": "group_level"})
            # Also compute mean norm_score and mean zscore per group
            grp_means = df_score.groupby(adata.obs[gcol].astype(str).values)[["norm_score", "zscore", "mc_pval"]].agg(
                mean_norm_score=("norm_score", "mean"),
                mean_zscore=("zscore", "mean"),
                prop_p_lt_005=("mc_pval", lambda x: np.mean(x < 0.05)),
            ).reset_index().rename(columns={"index": "group_level"})
            df_res["group_level"] = df_res["group_level"].astype(str)
            grp_means["group_level"] = grp_means["group_level"].astype(str)
            merged = pd.merge(df_res, grp_means, on="group_level", how="left")
            merged["atlas"] = atlas_key
            merged["trait"] = trait
            merged["group_col"] = gcol
            merged["n_genes_matched"] = len(g_list)
            all_group_stats.append(merged)

    if all_group_stats:
        df_out = pd.concat(all_group_stats, ignore_index=True)
        out_file = RESULTS_DIR / f"scdrs_group_stats_{atlas_key}.tsv"
        df_out.to_csv(out_file, sep="\t", index=False)
        print(f"  Saved scDRS group stats -> {out_file}")

    # Compute individual gene expression specificity across primary cell_type column
    primary_col = group_cols[0]
    core_genes = gs_dict["PSD_GWAS_Locus_Prioritized"][0] + ["WNT9A", "MON1B", "IFNA6", "MCEMP1", "CHST5", "RYR2"]
    core_genes = list(dict.fromkeys([g for g in core_genes if g in adata.var_names]))
    if core_genes:
        sub_x = adata[:, core_genes].X
        if hasattr(sub_x, "toarray"):
            sub_x = sub_x.toarray()
        df_expr = pd.DataFrame(sub_x, index=adata.obs[primary_col].astype(str).values, columns=core_genes)
        mean_expr = df_expr.groupby(level=0).mean()
        frac_expr = (df_expr > 0).groupby(level=0).mean()
        # Z-score across cell types per gene
        z_expr = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0) + 1e-8)

        records = []
        for ctype in mean_expr.index:
            n_cells = int((adata.obs[primary_col].astype(str) == ctype).sum())
            for g in core_genes:
                records.append({
                    "atlas": atlas_key,
                    "group_col": primary_col,
                    "cell_type": ctype,
                    "n_cells": n_cells,
                    "gene": g,
                    "mean_log1p_expr": round(float(mean_expr.loc[ctype, g]), 4),
                    "pct_expressing": round(100.0 * float(frac_expr.loc[ctype, g]), 2),
                    "specificity_z": round(float(z_expr.loc[ctype, g]), 4),
                })
        df_spec = pd.DataFrame(records)
        spec_file = RESULTS_DIR / f"gene_specificity_{atlas_key}.tsv"
        df_spec.to_csv(spec_file, sep="\t", index=False)
        print(f"  Saved gene specificity table -> {spec_file}")


if __name__ == "__main__":
    gs_dict = load_gs_file(DATA_DIR / "psd_scdrs_genesets.gs")
    for key in ["skin_fibroblast_subclusters", "skin_keratinocyte_subclusters", "skin_12datasets_integrated", "tabula_sapiens_skin", "scalp_trunk_skin_68k"]:
        url, desc = ATLAS_URLS[key]
        run_scdrs_on_atlas(key, url, desc, gs_dict)
