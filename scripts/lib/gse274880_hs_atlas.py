#!/usr/bin/env python3
"""
Helper module (formerly step 9). Hidradenitis Suppurativa (HS) vs. Pilonidal Sinus Disease (PSD) Single-Cell & Genetic Dissection
1. Download GEO GSE274880 (52 MB tar of 8 10x filtered_feature_bc_matrix.h5 files:
   2 Healthy Whole Skin, 1 Gluteal HS Hurley II, 5 Axillary HS Hurley II/III).
2. Annotate skin/HS cell types (including tunnel keratinocytes, papillary/Wnt+ fibroblasts,
   eccrine sweat gland, tertiary lymphoid B/plasma cells, mast cells, myeloid phagocytes).
3. Run scDRS (n_ctrl=5,000) across cell types and Disease Condition (Healthy vs. Gluteal HS vs. Axillary HS)
   and test differential expression of PSD risk genes in Gluteal/Axillary HS vs. Healthy Skin.
4. Query Open Targets for Hidradenitis Suppurativa (EFO_1000711) to quantify shared vs. divergent
   genetic architecture between PSD and HS.

Current role: only load_and_annotate_gse (GSE274880 atlas loading and annotation) is used, imported by step 03b.
Its own main run (early gene sets, Open Targets overlap) is superseded and its outputs were removed.
"""
import os
for _env_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_env_var] = "2"
import json
import tarfile
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
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parents[2]   # scripts/lib/ -> repo root
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = BASE_DIR / "figures"

SAMPLE_META = {
    "GSM8459837": ("Healthy_1", "Healthy", "Normal_Skin", "Healthy"),
    "GSM8459838": ("Healthy_2", "Healthy", "Normal_Skin", "Healthy"),
    "GSM8459839": ("HS_Axilla_1", "HS_Hurley_III", "Axilla", "HS_Axilla"),
    "GSM8459840": ("HS_Axilla_2", "HS_Hurley_III", "Axilla", "HS_Axilla"),
    "GSM8459841": ("HS_Axilla_3", "HS_Hurley_III", "Axilla", "HS_Axilla"),
    "GSM8459842": ("HS_Gluteal_4", "HS_Hurley_II", "Gluteal", "HS_Gluteal"),
    "GSM8459843": ("HS_Axilla_5", "HS_Hurley_II", "Axilla", "HS_Axilla"),
    "GSM8459844": ("HS_Axilla_6", "HS_Hurley_III", "Axilla", "HS_Axilla"),
}

CANONICAL_MARKERS = {
    "Dermal Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "Basal/Prolif Keratinocyte": ["KRT5", "KRT14", "TP63", "MKI67", "TOP2A"],
    "Diff/Tunnel Keratinocyte": ["KRT1", "KRT10", "KRT6A", "KRT16", "S100A8", "S100A9", "IVL"],
    "Eccrine/Sweat Gland": ["AQP5", "DCD", "SCGB2A2", "PIP", "KRT7", "CHRM3"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "ACKR1"],
    "Smooth Muscle / Pericyte": ["ACTA2", "TAGLN", "RGS5", "MYH11", "MCAM"],
    "T / NK Cell": ["CD3D", "CD3E", "TRAC", "NKG7", "GZMB"],
    "Plasma / B Cell (TLS)": ["MZB1", "JCHAIN", "IGHG1", "CD79A", "MS4A1"],
    "Macrophage / Myeloid DC": ["LYZ", "CD68", "C1QA", "CD14", "HLA-DRA", "FCER1A"],
    "Mast Cell": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
    "Melanocyte / Schwann": ["PMEL", "MLANA", "DCT", "TYRP1", "SOX10"],
}


def load_and_annotate_gse274880() -> ad.AnnData:
    gse_dir = DATA_DIR / "GSE274880"
    gse_dir.mkdir(exist_ok=True)
    tar_path = gse_dir / "GSE274880_RAW.tar"
    if not tar_path.exists():
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274880/suppl/GSE274880_RAW.tar"
        print(f"[1/4] Downloading GSE274880 (52 MB) -> {tar_path}...")
        urllib.request.urlretrieve(url, tar_path)
    if not list(gse_dir.glob("*.h5")):
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=gse_dir)

    adatas = []
    for h5_file in sorted(gse_dir.glob("*.h5")):
        gsm = h5_file.name.split("_")[0]
        if gsm not in SAMPLE_META:
            continue
        s_label, stage, site, group_cond = SAMPLE_META[gsm]
        a = sc.read_10x_h5(h5_file)
        a.var_names_make_unique()
        sc.pp.filter_cells(a, min_genes=300)
        a.obs["gsm"] = gsm
        a.obs["sample_label"] = s_label
        a.obs["disease_stage"] = stage
        a.obs["anatomical_site"] = site
        a.obs["condition_group"] = group_cond
        adatas.append(a)

    adata = ad.concat(adatas, join="outer", index_unique="-")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs["n_genes_by_counts"] < 7000) & (adata.obs["pct_counts_mt"] < 20.0)].copy()
    print(f"  Combined QC-filtered GSE274880 atlas: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(adata.obs["condition_group"].value_counts())

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, subset=False)
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden")

    # Score canonical marker sets (Z-standardized across cells) to assign biological cell types to Leiden clusters
    score_cols = []
    for ctype, markers in CANONICAL_MARKERS.items():
        present = [m for m in markers if m in adata.var_names]
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ctype}")
            vals = adata.obs[f"score_{ctype}"].values
            adata.obs[f"score_{ctype}"] = (vals - np.mean(vals)) / (np.std(vals) + 1e-8)
            score_cols.append((ctype, f"score_{ctype}"))

    cluster_map = {}
    for cl, grp in adata.obs.groupby("leiden"):
        best_ct = max(score_cols, key=lambda x: grp[x[1]].mean())[0]
        cluster_map[cl] = best_ct

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_map).astype(str)
    adata.obs["cell_type_condition"] = adata.obs["cell_type"] + " | " + adata.obs["condition_group"].astype(str)
    print("  Annotated cell type counts in GSE274880:")
    print(adata.obs["cell_type"].value_counts())
    return adata


def run_scdrs_on_hs(adata: ad.AnnData):
    print("\n[2/4] Running scDRS (n_ctrl=1,000) on GSE274880 (Healthy vs Gluteal HS vs Axillary HS)...")
    scdrs.preprocess(adata, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)

    gs_dict = {}
    with open(DATA_DIR / "psd_scdrs_genesets.gs") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            t, gstr = line.strip().split("\t")
            pairs = [x.split(":") for x in gstr.split(",") if ":" in x]
            gs_dict[t] = ([p[0] for p in pairs], [float(p[1]) for p in pairs])

    all_rows = []
    for trait, (g_list, w_list) in gs_dict.items():
        if trait == "PSD_Polygenic_Top1000" or trait.startswith("PSD_Module_"):
            continue
        present = [(g, w) for g, w in zip(g_list, w_list) if g in adata.var_names]
        if len(present) < 3:
            continue
        g_p = [x[0] for x in present]
        w_p = [x[1] for x in present]
        print(f"  Scoring {trait} ({len(g_p)} genes)...")
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

        for gcol in ["cell_type", "condition_group", "cell_type_condition"]:
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
            ).reset_index().rename(columns={"index": "group_level"})
            df_res["group_level"] = df_res["group_level"].astype(str)
            grp_means["group_level"] = grp_means["group_level"].astype(str)
            merged = pd.merge(df_res, grp_means, on="group_level", how="left")
            merged["atlas"] = "GSE274880_HS_and_Healthy_Skin"
            merged["trait"] = trait
            merged["group_col"] = gcol
            merged["n_genes_matched"] = len(g_p)
            merged["pct_fdr_0.1"] = (100.0 * merged["n_fdr_0.1"] / merged["n_cell"]).round(2)
            all_rows.append(merged)

    df_out = pd.concat(all_rows, ignore_index=True)
    out_file = RESULTS_DIR / "scdrs_group_stats_GSE274880_HS.tsv"
    df_out.sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False]).to_csv(out_file, sep="\t", index=False)
    print(f"  Saved GSE274880 HS scDRS summary -> {out_file}")
    top = df_out[(df_out["assoc_mcp"] < 0.05) & (df_out["n_cell"] >= 30)].sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False])
    cols = ["trait", "group_col", "group_level", "n_cell", "assoc_mcp", "assoc_mcz", "hetero_mcp", "prop_p_lt_005", "pct_fdr_0.1"]
    print(top[cols].head(35).to_string(index=False))

    # Plot heatmap of cell_type_condition scDRS assoc_mcz for Dermal Fibroblast, Keratinocytes, Eccrine, Plasma/B, Myeloid
    sub_ct = df_out[df_out["group_col"] == "cell_type"].copy()
    piv_z = sub_ct.pivot(index="group_level", columns="trait", values="assoc_mcz").fillna(0.0)
    piv_p = sub_ct.pivot(index="group_level", columns="trait", values="assoc_mcp").fillna(1.0)
    annot = piv_z.round(2).astype(str)
    for r in piv_z.index:
        for c in piv_z.columns:
            p = piv_p.loc[r, c]
            star = "**" if p < 0.01 else ("*" if p < 0.05 else ("+" if p < 0.10 else ""))
            annot.loc[r, c] = f"{piv_z.loc[r, c]:.2f}{star}"
    fig, ax = plt.subplots(figsize=(12, 6.5))
    sns.heatmap(piv_z, cmap="RdBu_r", center=0, vmin=-3, vmax=4.5, annot=annot, fmt="", linewidths=0.5, ax=ax)
    ax.set_title("scDRS Enrichment of PSD Genetic Sets Across Hidradenitis Suppurativa & Healthy Skin (GSE274880)", fontweight="bold")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "scdrs_heatmap_GSE274880_HS.png", dpi=300)
    plt.close(fig)

    # Compare expression of PSD core genes in Fibroblasts & Keratinocytes across Healthy vs Gluteal HS vs Axillary HS
    print("\n[3/4] Comparing PSD Risk Gene Expression across Healthy Skin vs. Gluteal HS vs. Axillary HS...")
    core_genes = [g for g in ["TWIST1", "HDAC9", "DKK1", "TBX15", "WARS2", "PDGFA", "WNT9A", "PRKAR1B", "NEK9", "TMED10", "TGM3", "PADI3", "TCHH", "CHRM3", "TTK", "MON1B", "MCEMP1", "RNFT2", "MBL2"] if g in adata.var_names]
    sub_x = adata[:, core_genes].X
    if hasattr(sub_x, "toarray"):
        sub_x = sub_x.toarray()
    df_e = pd.DataFrame(sub_x, columns=core_genes)
    df_e["cell_type"] = adata.obs["cell_type"].values
    df_e["condition_group"] = adata.obs["condition_group"].values
    grp_e = df_e.groupby(["cell_type", "condition_group"])[core_genes].mean().reset_index()
    grp_e.to_csv(RESULTS_DIR / "gse274880_hs_vs_healthy_psd_gene_expression.tsv", sep="\t", index=False)
    for ct in ["Dermal Fibroblast", "Basal/Prolif Keratinocyte", "Diff/Tunnel Keratinocyte", "Eccrine/Sweat Gland", "Plasma / B Cell (TLS)"]:
        sub_c = grp_e[grp_e["cell_type"] == ct]
        if not sub_c.empty:
            print(f"\n  --- {ct} (Mean log1p expression: Healthy vs HS_Gluteal vs HS_Axilla) ---")
            print(sub_c.set_index("condition_group")[["TWIST1", "HDAC9", "DKK1", "TBX15", "WARS2", "PDGFA", "WNT9A", "NEK9", "TMED10", "TGM3", "CHRM3", "TTK", "MON1B"]].round(4).to_string())


def compare_psd_vs_hs_opentargets():
    print("\n[4/4] Comparing Open Targets Genetic Architecture: Pilonidal Sinus/Abscess vs. Hidradenitis Suppurativa (MONDO_0006559)...")
    query = """
    query($efoId: String!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: 0, size: 500}) {
          count
          rows {
            score
            target {
              id
              approvedSymbol
              approvedName
            }
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """
    req = urllib.request.Request(
        "https://api.platform.opentargets.org/api/v4/graphql",
        data=json.dumps({"query": query, "variables": {"efoId": "MONDO_0006559"}}).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"}
    )
    res = json.loads(urllib.request.urlopen(req, timeout=30).read().decode("utf-8"))
    rows = (res.get("data", {}).get("disease") or {}).get("associatedTargets", {}).get("rows", [])
    hs_records = []
    for r in rows:
        dtypes = {x["id"]: x["score"] for x in r.get("datatypeScores", [])}
        hs_records.append({
            "gene_symbol": r["target"]["approvedSymbol"],
            "ensembl_id": r["target"]["id"],
            "gene_name": r["target"]["approvedName"],
            "hs_overall_score": r["score"],
            "hs_genetic_score": dtypes.get("genetic_association", 0.0),
            "hs_somatic_or_rna_score": dtypes.get("rna_expression", 0.0),
            "hs_drug_score": dtypes.get("known_drug", 0.0),
            "hs_literature_score": dtypes.get("literature", 0.0),
        })
    df_hs = pd.DataFrame(hs_records)
    df_hs.to_csv(DATA_DIR / "opentargets_hidradenitis_targets.tsv", sep="\t", index=False)

    # Compare against PSD targets (OpenTargets + Curated + Exome Top 100)
    df_psd_ot = pd.read_csv(DATA_DIR / "opentargets_pilonidal_targets.tsv", sep="\t")
    df_psd_cur = pd.read_csv(RESULTS_DIR / "curated_psd_risk_genes.tsv", sep="\t")
    df_psd_poly = pd.read_csv(RESULTS_DIR / "psd_genomewide_polygenic_gene_scores.tsv", sep="\t").head(200)

    psd_genes = set(df_psd_ot["gene_symbol"]).union(set(df_psd_cur["gene_symbol"])).union(set(df_psd_poly["gene_symbol"]))
    overlap = df_hs[df_hs["gene_symbol"].isin(psd_genes)].copy()
    overlap = pd.merge(
        overlap,
        df_psd_poly[["gene_symbol", "z_polygenic_combined", "z_exome", "z_mvp_gwas_100kb", "z_opentargets"]],
        on="gene_symbol",
        how="left"
    ).sort_values(["hs_genetic_score", "hs_overall_score"], ascending=False)
    out_ov = RESULTS_DIR / "psd_vs_hs_shared_targets.tsv"
    overlap.to_csv(out_ov, sep="\t", index=False)
    print(f"  Saved {len(overlap)} shared PSD-HS OpenTargets genes -> {out_ov}")
    print(overlap.head(25).to_string(index=False))


def clean_and_replot_cached_hs():
    out_file = RESULTS_DIR / "scdrs_group_stats_GSE274880_HS.tsv"
    df_out = pd.read_csv(out_file, sep="\t")
    df_out = df_out[~df_out["trait"].astype(str).str.startswith("PSD_Module_")].copy()
    df_out.sort_values(["assoc_mcp", "assoc_mcz"], ascending=[True, False]).to_csv(out_file, sep="\t", index=False)
    sub_ct = df_out[df_out["group_col"] == "cell_type"].copy()
    piv_z = sub_ct.pivot(index="group_level", columns="trait", values="assoc_mcz").fillna(0.0)
    piv_p = sub_ct.pivot(index="group_level", columns="trait", values="assoc_mcp").fillna(1.0)
    annot = piv_z.round(2).astype(str)
    for r in piv_z.index:
        for c in piv_z.columns:
            p = piv_p.loc[r, c]
            star = "**" if p < 0.01 else ("*" if p < 0.05 else ("+" if p < 0.10 else ""))
            annot.loc[r, c] = f"{piv_z.loc[r, c]:.2f}{star}"
    fig, ax = plt.subplots(figsize=(11, 6.2))
    sns.heatmap(piv_z, cmap="RdBu_r", center=0, vmin=-3, vmax=4.5, annot=annot, fmt="", linewidths=0.5, ax=ax)
    ax.set_title("scDRS Enrichment of Unbiased PSD Genetic Sets Across Hidradenitis Suppurativa (GSE274880)", fontweight="bold")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "scdrs_heatmap_GSE274880_HS.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    if (RESULTS_DIR / "scdrs_group_stats_GSE274880_HS.tsv").exists() and (RESULTS_DIR / "gse274880_hs_vs_healthy_psd_gene_expression.tsv").exists():
        clean_and_replot_cached_hs()
    else:
        adata_hs = load_and_annotate_gse274880()
        run_scdrs_on_hs(adata_hs)
    compare_psd_vs_hs_opentargets()

