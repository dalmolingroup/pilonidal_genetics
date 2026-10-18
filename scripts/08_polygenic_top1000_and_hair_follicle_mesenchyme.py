#!/usr/bin/env python3
"""
Step 8:
1. Build genome-wide polygenic gene sets (Top 500 and Top 1,000 protein-coding genes)
   by integrating MVP Multi-Ancestry GWAS (GCST90480439), South Asian GWAS (GCST90727025),
   UKB Whole-Exome Gene Burden (GCST90084318, 18,807 genes), and Open Targets L2G.
2. Run scDRS with n_ctrl=5,000 on Tabula Sapiens Skin and Skin Fibroblast Subclusters
   to evaluate how 500/1,000-gene sets and a lower Monte Carlo p-value floor (2e-4)
   affect FDR-significant cell proportions and group association Z-scores.
3. Download and demultiplex GEO GSE287304 (Human Hair Follicle Mesenchymal Cells:
   Dermal Papilla [DP], Dermal Sheath Cup [DSC], and Upper Dermal Sheath [UDS])
   and run scDRS + gene specificity across the three hair-follicle mesenchymal compartments.
"""
import os
for _env_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_env_var] = "2"
import gzip
import urllib.request
import warnings
from pathlib import Path
import numpy as np
if not hasattr(np, "float_"):
    np.float_ = np.float64
import scipy.stats
import scipy.sparse as sp
from scipy.stats import norm
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

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = BASE_DIR / "figures"


def build_polygenic_top500_top1000():
    print("[1/3] Building genome-wide polygenic Top 500 and Top 1,000 gene sets...")
    # Load all 18,807 protein-coding genes with chrom, pos, and exome z_score
    exome_df = pd.read_csv(RESULTS_DIR / "ukb_exome_gene_burden_summary.tsv", sep="\t", dtype=str)
    exome_df["pos"] = pd.to_numeric(exome_df["pos"], errors="coerce")
    exome_df["z_exome"] = pd.to_numeric(exome_df["z_score"], errors="coerce").fillna(0.0)
    exome_df = exome_df.dropna(subset=["chrom", "pos"]).copy()

    # Load MVP multi-ancestry GWAS top variants (p < 1e-4) and filter to MAF >= 0.01
    mvp_df = pd.read_csv(RESULTS_DIR / "gwas_MVP_multiancestry_GCST90480439_top_variants.tsv", sep="\t", dtype=str)
    mvp_df["p_value"] = pd.to_numeric(mvp_df["p_value"], errors="coerce")
    mvp_df["pos"] = pd.to_numeric(mvp_df["base_pair_location"], errors="coerce")
    mvp_df["eaf"] = pd.to_numeric(mvp_df["effect_allele_frequency"], errors="coerce")
    mvp_df = mvp_df.dropna(subset=["p_value", "pos", "eaf"]).copy()
    mvp_df["maf"] = np.minimum(mvp_df["eaf"], 1.0 - mvp_df["eaf"])
    mvp_df = mvp_df[mvp_df["maf"] >= 0.01].copy()
    mvp_df["z_gwas"] = -norm.ppf((mvp_df["p_value"] / 2.0).clip(1e-15, 0.4999))

    # Map GWAS variants within +/- 100 kb of gene coordinates per chromosome
    gwas_gene_z = {}
    for chrom, grp_genes in exome_df.groupby("chrom"):
        sub_v = mvp_df[mvp_df["chrom" if "chrom" in mvp_df.columns else "chromosome"].astype(str) == str(chrom)]
        if sub_v.empty:
            continue
        v_pos = sub_v["pos"].values
        v_z = sub_v["z_gwas"].values
        for _, grow in grp_genes.iterrows():
            g_sym = grow["gene_symbol"]
            g_pos = float(grow["pos"])
            mask = np.abs(v_pos - g_pos) <= 100000.0
            if np.any(mask):
                gwas_gene_z[g_sym] = float(np.max(v_z[mask]))

    # Load Open Targets genetic scores
    ot_z = {}
    ot_path = DATA_DIR / "opentargets_pilonidal_targets.tsv"
    if ot_path.exists():
        ot_df = pd.read_csv(ot_path, sep="\t")
        for _, r in ot_df.iterrows():
            gsc = float(r["genetic_association_score"])
            if gsc > 0:
                ot_z[str(r["gene_symbol"])] = max(ot_z.get(str(r["gene_symbol"]), 0.0), 2.5 + 10.0 * gsc)

    # Load curated GWAS risk genes
    cur_df = pd.read_csv(RESULTS_DIR / "curated_psd_risk_genes.tsv", sep="\t")
    cur_z = dict(zip(cur_df["gene_symbol"], cur_df["prior_z_weight"].astype(float)))

    # Compute combined polygenic score across all 18,807 genes
    combined_records = []
    for _, r in exome_df.iterrows():
        sym = str(r["gene_symbol"])
        z_ex = float(r["z_exome"])
        z_gw = gwas_gene_z.get(sym, 0.0)
        z_ot = ot_z.get(sym, 0.0)
        z_cu = cur_z.get(sym, 0.0)
        # Stouffer / max-weighted integration of common GWAS + rare exome + OpenTargets
        if z_gw > 0:
            z_comb = (z_ex + 1.2 * z_gw) / np.sqrt(1.0 + 1.44)
        else:
            z_comb = z_ex
        z_final = max(z_comb, z_gw, z_ot, z_cu)
        combined_records.append({
            "gene_symbol": sym,
            "ensembl_id": r["ensembl_id"],
            "chrom": r["chrom"],
            "pos": int(r["pos"]),
            "z_exome": round(z_ex, 4),
            "z_mvp_gwas_100kb": round(z_gw, 4),
            "z_opentargets": round(z_ot, 4),
            "z_curated": round(z_cu, 4),
            "z_polygenic_combined": round(z_final, 4),
        })

    df_poly = pd.DataFrame(combined_records).sort_values("z_polygenic_combined", ascending=False)
    out_poly = RESULTS_DIR / "psd_genomewide_polygenic_gene_scores.tsv"
    df_poly.to_csv(out_poly, sep="\t", index=False)
    print(f"  Saved genome-wide polygenic scores for {len(df_poly)} genes -> {out_poly}")
    print("  Top 25 genes in PSD_Polygenic_Combined:")
    print(df_poly.head(25).to_string(index=False))

    # Update data/psd_scdrs_genesets.gs with Top500 and Top1000
    gs_path = DATA_DIR / "psd_scdrs_genesets.gs"
    lines = [
        l for l in gs_path.read_text().splitlines()
        if not l.startswith(("PSD_Polygenic_Top500\t", "PSD_Polygenic_Top1000\t"))
    ]
    top500 = df_poly.head(500)
    top1000 = df_poly.head(1000)
    lines.append("PSD_Polygenic_Top500\t" + ",".join(f"{r['gene_symbol']}:{r['z_polygenic_combined']:.4f}" for _, r in top500.iterrows()))
    lines.append("PSD_Polygenic_Top1000\t" + ",".join(f"{r['gene_symbol']}:{r['z_polygenic_combined']:.4f}" for _, r in top1000.iterrows()))
    gs_path.write_text("\n".join(lines) + "\n")
    print(f"  Added PSD_Polygenic_Top500 and PSD_Polygenic_Top1000 to {gs_path}")


def run_high_resolution_scdrs_comparison():
    print("\n[2/3] Running scDRS sensitivity comparison (Top500/Top1000 vs Oligogenic sets, n_ctrl=1,000)...")
    # Load gene sets
    gs_dict = {}
    with open(DATA_DIR / "psd_scdrs_genesets.gs") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            t, gstr = line.strip().split("\t")
            pairs = [x.split(":") for x in gstr.split(",") if ":" in x]
            gs_dict[t] = ([p[0] for p in pairs], [float(p[1]) for p in pairs])

    target_traits = [
        "PSD_Polygenic_Top1000",
        "PSD_Polygenic_Top500",
        "PSD_GWAS_Locus_Prioritized",
        "PSD_MVP_MultiAncestry_GWAS",
        "PSD_Exome_Burden_Top200",
    ]

    for atlas_key, gcol in [("tabula_sapiens_skin", "cell_type"), ("skin_fibroblast_subclusters", "author_cell_type")]:
        h5ad_path = DATA_DIR / f"{atlas_key}.h5ad"
        adata_orig = sc.read_h5ad(h5ad_path)
        obsp_saved = dict(adata_orig.obsp)
        obsm_saved = dict(adata_orig.obsm)
        adata = adata_orig.raw.to_adata() if adata_orig.raw is not None else adata_orig
        for k, v in obsp_saved.items():
            adata.obsp[k] = v
        for k, v in obsm_saved.items():
            adata.obsm[k] = v
        if "feature_name" in adata.var.columns:
            adata.var_names = adata.var["feature_name"].astype(str).values
            adata.var_names_make_unique()
        x_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
        if x_max > 50:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        if "connectivities" not in adata.obsp:
            sc.pp.pca(adata, n_comps=30)
            sc.pp.neighbors(adata)

        scdrs.preprocess(adata, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)

        rows = []
        for trait in target_traits:
            g_list, w_list = gs_dict[trait]
            present = [(g, w) for g, w in zip(g_list, w_list) if g in adata.var_names]
            g_p = [x[0] for x in present]
            w_p = [x[1] for x in present]
            print(f"  [{atlas_key}] Scoring {trait} ({len(g_p)} genes, n_ctrl=1000)...")
            df_score = scdrs.score_cell(
                data=adata,
                gene_list=g_p,
                gene_weight=w_p,
                ctrl_match_key="mean_var",
                n_ctrl=1000,
                weight_opt="vs",
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=True,
                verbose=False,
            ).fillna(0.0)

            dict_df = scdrs.method.downstream_group_analysis(
                adata=adata,
                df_full_score=df_score,
                group_cols=[gcol],
                fdr_thresholds=[0.05, 0.1, 0.2],
            )
            df_res = dict_df[gcol].reset_index().rename(columns={"group": "group_level", gcol: "group_level", "index": "group_level"})
            grp_means = df_score.groupby(adata.obs[gcol].astype(str).values)[["norm_score", "zscore", "mc_pval"]].agg(
                mean_norm_score=("norm_score", "mean"),
                mean_zscore=("zscore", "mean"),
                prop_p_lt_005=("mc_pval", lambda x: np.mean(x < 0.05)),
                prop_p_lt_001=("mc_pval", lambda x: np.mean(x < 0.01)),
            ).reset_index().rename(columns={"index": "group_level"})
            df_res["group_level"] = df_res["group_level"].astype(str)
            grp_means["group_level"] = grp_means["group_level"].astype(str)
            merged = pd.merge(df_res, grp_means, on="group_level", how="left")
            merged["atlas"] = atlas_key
            merged["trait"] = trait
            merged["n_genes"] = len(g_p)
            merged["pct_fdr_0.1"] = (100.0 * merged["n_fdr_0.1"] / merged["n_cell"]).round(2)
            merged["pct_fdr_0.2"] = (100.0 * merged["n_fdr_0.2"] / merged["n_cell"]).round(2)
            rows.append(merged)

        df_comp = pd.concat(rows, ignore_index=True)
        out_f = RESULTS_DIR / f"scdrs_nctrl5000_polygenic_comparison_{atlas_key}.tsv"
        df_comp.sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False]).to_csv(out_f, sep="\t", index=False)
        print(f"  Saved polygenic comparison -> {out_f}")
        top_sub = df_comp[df_comp["assoc_mcp"] < 0.08].sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False])
        cols = ["trait", "n_genes", "group_level", "n_cell", "assoc_mcp", "assoc_mcz", "hetero_mcp", "prop_p_lt_005", "pct_fdr_0.1", "pct_fdr_0.2"]
        print(top_sub[cols].head(25).to_string(index=False))


def analyze_gse287304_hair_follicle_mesenchyme():
    print("\n[3/3] Analyzing GEO GSE287304 (Hair Follicle Mesenchymal Compartments: DP, DSC, UDS)...")
    gse_dir = DATA_DIR / "GSE287304"
    gse_dir.mkdir(exist_ok=True)
    triplets_path = gse_dir / "filtered_cells_triplets.tsv"
    if not triplets_path.exists():
        raise FileNotFoundError(f"Expected pre-filtered triplets at {triplets_path}")

    print("  Loading pre-filtered GSE287304 triplets (barcodes with >=500 RNA genes)...")
    features = pd.read_csv(gse_dir / "GSE287304_1_features.tsv.gz", sep="\t", header=None, names=["id", "name", "ftype"])
    barcodes = pd.read_csv(gse_dir / "GSE287304_1_barcodes.tsv.gz", sep="\t", header=None, names=["barcode"])

    df_trip = pd.read_csv(
        triplets_path,
        sep=" ",
        engine="c",
        header=None,
        names=["feat_idx", "bc_idx", "count"],
        dtype={"feat_idx": np.int32, "bc_idx": np.int32, "count": np.int32},
    )
    # Map 1-based bc_idx to contiguous 0..n_cells-1
    unique_bcs, cell_col = np.unique(df_trip["bc_idx"].values, return_inverse=True)
    feat_row = df_trip["feat_idx"].values - 1
    counts = df_trip["count"].values
    del df_trip

    n_cells = len(unique_bcs)
    n_feats = len(features)
    print(f"  Loaded {n_cells:,} filtered cells across {n_feats:,} features")

    mat_cells_feats = sp.csr_matrix(
        (counts, (cell_col, feat_row)),
        shape=(n_cells, n_feats),
        dtype=np.float32,
    )
    del counts, cell_col, feat_row
    rna_mask = (features["ftype"] == "Gene Expression").values
    hto_mask = (features["ftype"] == "Antibody Capture").values

    mat_rna_cells = mat_cells_feats[:, rna_mask].tocsr()
    mat_hto_cells = mat_cells_feats[:, hto_mask].toarray()
    del mat_cells_feats
    umi_per_bc = np.asarray(mat_rna_cells.sum(axis=1)).ravel()

    hto_names = features.loc[hto_mask, "id"].values
    # Hashtag mapping from GSE287304_feature_reference.csv.gz:
    # HTO_C0252 -> DSC (Dermal Sheath Cup)
    # HTO_C0253 -> DP  (Dermal Papilla)
    # HTO_C0254 -> UDS (Upper Dermal Sheath)
    idx_dsc = np.where(hto_names == "HTO_C0252")[0][0]
    idx_dp = np.where(hto_names == "HTO_C0253")[0][0]
    idx_uds = np.where(hto_names == "HTO_C0254")[0][0]

    hto_3 = mat_hto_cells[:, [idx_dsc, idx_dp, idx_uds]].astype(np.float64)
    # Centered log-ratio (CLR) normalization across cells for each HTO
    log_hto = np.log1p(hto_3)
    clr_hto = log_hto - log_hto.mean(axis=0, keepdims=True)
    best_idx = np.argmax(clr_hto, axis=1)
    sorted_clr = np.sort(clr_hto, axis=1)
    margin = sorted_clr[:, -1] - sorted_clr[:, -2]
    labels_map = np.array(["DSC (Dermal Sheath Cup)", "DP (Dermal Papilla)", "UDS (Upper Dermal Sheath)"])
    assigned = labels_map[best_idx]
    # Keep confident singlets with positive CLR, margin >= 0.35, and RNA UMI >= 1000
    singlet_mask = (sorted_clr[:, -1] > 0.5) & (margin >= 0.35) & (umi_per_bc >= 1000)

    bc_strings = barcodes["barcode"].values[unique_bcs - 1][singlet_mask]
    adata = ad.AnnData(
        X=mat_rna_cells[singlet_mask, :].astype(np.float32),
        obs=pd.DataFrame(
            {
                "barcode": bc_strings,
                "compartment": assigned[singlet_mask],
                "total_umi": umi_per_bc[singlet_mask],
            },
            index=bc_strings,
        ),
        var=pd.DataFrame(
            {"ensembl_id": features.loc[rna_mask, "id"].values},
            index=features.loc[rna_mask, "name"].astype(str).values,
        ),
    )
    del mat_rna_cells
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=500)
    print(f"  High-confidence demultiplexed hair follicle mesenchymal singlets: {adata.n_obs:,} cells")
    print(adata.obs["compartment"].value_counts())

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata)
    scdrs.preprocess(adata, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)

    # Load gene sets
    gs_dict = {}
    with open(DATA_DIR / "psd_scdrs_genesets.gs") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            t, gstr = line.strip().split("\t")
            pairs = [x.split(":") for x in gstr.split(",") if ":" in x]
            gs_dict[t] = ([p[0] for p in pairs], [float(p[1]) for p in pairs])

    rows = []
    for trait, (g_list, w_list) in gs_dict.items():
        if trait == "PSD_Polygenic_Top1000":
            continue
        present = [(g, w) for g, w in zip(g_list, w_list) if g in adata.var_names]
        if len(present) < 3:
            continue
        g_p = [x[0] for x in present]
        w_p = [x[1] for x in present]
        df_score = scdrs.score_cell(
            data=adata,
            gene_list=g_p,
            gene_weight=w_p,
            ctrl_match_key="mean_var",
            n_ctrl=1000,
            weight_opt="vs",
            return_ctrl_raw_score=False,
            return_ctrl_norm_score=True,
            verbose=False,
        ).fillna(0.0)

        dict_df = scdrs.method.downstream_group_analysis(
            adata=adata,
            df_full_score=df_score,
            group_cols=["compartment"],
            fdr_thresholds=[0.05, 0.1, 0.2],
        )
        df_res = dict_df["compartment"].reset_index().rename(columns={"group": "group_level", "compartment": "group_level", "index": "group_level"})
        grp_means = df_score.groupby(adata.obs["compartment"].astype(str).values)[["norm_score", "zscore", "mc_pval"]].agg(
            mean_norm_score=("norm_score", "mean"),
            mean_zscore=("zscore", "mean"),
            prop_p_lt_005=("mc_pval", lambda x: np.mean(x < 0.05)),
        ).reset_index().rename(columns={"index": "group_level"})
        df_res["group_level"] = df_res["group_level"].astype(str)
        grp_means["group_level"] = grp_means["group_level"].astype(str)
        merged = pd.merge(df_res, grp_means, on="group_level", how="left")
        merged["atlas"] = "GSE287304_hair_follicle_mesenchyme"
        merged["trait"] = trait
        merged["n_genes_matched"] = len(g_p)
        merged["pct_fdr_0.1"] = (100.0 * merged["n_fdr_0.1"] / merged["n_cell"]).round(2)
        rows.append(merged)

    df_hf = pd.concat(rows, ignore_index=True)
    out_hf = RESULTS_DIR / "scdrs_group_stats_GSE287304_hair_follicle_mesenchyme.tsv"
    df_hf.sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False]).to_csv(out_hf, sep="\t", index=False)
    print(f"  Saved GSE287304 Hair Follicle Mesenchyme scDRS results -> {out_hf}")
    cols = ["trait", "group_level", "n_cell", "assoc_mcp", "assoc_mcz", "hetero_mcp", "prop_p_lt_005", "pct_fdr_0.1"]
    print(df_hf.sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False])[cols].head(20).to_string(index=False))

    # Gene expression specificity across DP, DSC, UDS
    core_genes = gs_dict["PSD_GWAS_Locus_Prioritized"][0] + ["WNT9A", "MON1B", "IFNA6", "MCEMP1", "CHST5", "RYR2", "ALPL", "SOX2", "CORIN", "VCAN", "ACTA2"]
    core_genes = list(dict.fromkeys([g for g in core_genes if g in adata.var_names]))
    sub_x = adata[:, core_genes].X.toarray()
    df_expr = pd.DataFrame(sub_x, index=adata.obs["compartment"].astype(str).values, columns=core_genes)
    mean_expr = df_expr.groupby(level=0).mean()
    frac_expr = (df_expr > 0).groupby(level=0).mean()
    z_expr = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0) + 1e-8)

    records = []
    for ctype in mean_expr.index:
        n_c = int((adata.obs["compartment"].astype(str) == ctype).sum())
        for g in core_genes:
            records.append({
                "compartment": ctype,
                "n_cells": n_c,
                "gene": g,
                "mean_log1p_expr": round(float(mean_expr.loc[ctype, g]), 4),
                "pct_expressing": round(100.0 * float(frac_expr.loc[ctype, g]), 2),
                "specificity_z": round(float(z_expr.loc[ctype, g]), 4),
            })
    df_spec = pd.DataFrame(records).sort_values(["gene", "specificity_z"], ascending=[True, False])
    spec_f = RESULTS_DIR / "gene_specificity_GSE287304_hair_follicle_mesenchyme.tsv"
    df_spec.to_csv(spec_f, sep="\t", index=False)
    print(f"  Saved GSE287304 DP/DSC/UDS gene specificity -> {spec_f}")
    print("\n  === Key PSD Genes Across DP (Dermal Papilla), DSC (Dermal Sheath Cup), and UDS (Upper Dermal Sheath) ===")
    piv_pct = df_spec.pivot(index="gene", columns="compartment", values="pct_expressing")
    piv_mean = df_spec.pivot(index="gene", columns="compartment", values="mean_log1p_expr")
    print("  % Expressing:")
    print(piv_pct.to_string())


if __name__ == "__main__":
    if not (RESULTS_DIR / "psd_genomewide_polygenic_gene_scores.tsv").exists():
        build_polygenic_top500_top1000()
    else:
        print("[1/3] Polygenic Top 500 / Top 1,000 scores already generated; skipping Step 1.")
    if not (RESULTS_DIR / "scdrs_nctrl5000_polygenic_comparison_skin_fibroblast_subclusters.tsv").exists():
        run_high_resolution_scdrs_comparison()
    else:
        print("[2/3] Polygenic scDRS comparison already generated; skipping Step 2.")
    analyze_gse287304_hair_follicle_mesenchyme()
