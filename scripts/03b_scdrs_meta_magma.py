#!/usr/bin/env python3
"""
Step 03b: scDRS (n_ctrl=1,000; z-score-weighted top-1,000 MAGMA genes) for PSD meta-analysis gene sets and comparator traits
across the 5 CELLxGENE skin atlases (atlas helpers in scripts/lib/scdrs_atlases.py), the GSE274880 HS/healthy skin atlas (scripts/lib/gse274880_hs_atlas.py) and the GSE287304
hair-follicle mesenchyme atlas (HTO demultiplexing rebuilt from the public GEO matrix).
Gene sets: data/psd_meta_magma_scdrs.gs (from 03a). Outputs: results/scdrs_magma_group_stats_<atlas>.tsv, results/scdrs_magma_all_atlases.tsv
Environment overrides: SCDRS_GS (gene-set file, e.g. the MHC-excluded one from 03c), SCDRS_TAG (output prefix, default scdrs_magma),
SCDRS_ATLAS (score only this atlas key, for SLURM arrays), SCDRS_MERGE=1 (only merge per-atlas tables and compute FDR; no scoring).
"""
import importlib.util, os, sys, urllib.request, gzip
from pathlib import Path
import numpy as np
BASE = Path(__file__).resolve().parent.parent
TAG = os.environ.get("SCDRS_TAG", "scdrs_magma")
ATLAS_KEYS = ["tabula_sapiens_skin", "skin_12datasets_integrated", "skin_keratinocyte_subclusters", "skin_fibroblast_subclusters",
              "scalp_trunk_skin_68k", "GSE274880_HS", "GSE287304_hair_follicle_mesenchyme"]


def merge_atlases():
    # BH-FDR within each trait across all atlas x group tests, from assoc_mcp (floored at 1/(n_ctrl+1)) and from assoc_mcz
    import pandas as pd
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
    res = BASE / "results"
    allr = pd.concat([pd.read_csv(res / f"{TAG}_group_stats_{k}.tsv", sep="\t") for k in ATLAS_KEYS], ignore_index=True)
    allr["assoc_fdr_within_trait"] = allr.groupby("trait")["assoc_mcp"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    allr["assoc_p_mcz"] = norm.sf(allr["assoc_mcz"])
    allr["assoc_fdr_mcz_within_trait"] = allr.groupby("trait")["assoc_p_mcz"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    allr.to_csv(res / f"{TAG}_all_atlases.tsv", sep="\t", index=False)
    print(allr.sort_values("assoc_p_mcz").groupby("trait").head(8)[["trait", "atlas", "group_level", "n_cell", "assoc_mcp", "assoc_mcz",
          "assoc_fdr_mcz_within_trait", "hetero_mcp", "n_fdr_0.1"]].to_string(index=False))


if os.environ.get("SCDRS_MERGE") == "1":
    merge_atlases(); print("SCDRS_MAGMA_DONE"); sys.exit(0)


def load_module(name, fname):
    spec = importlib.util.spec_from_file_location(name, BASE / "scripts" / fname)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
s03 = load_module("s03", "lib/scdrs_atlases.py")        # also applies numpy/scipy compatibility patches
s09 = load_module("s09", "lib/gse274880_hs_atlas.py")
import pandas as pd, scanpy as sc, anndata as ad, scdrs, scipy.sparse as sp, scipy.io
DATA, RES = BASE / "data", BASE / "results"
gs = s03.load_gs_file(Path(os.environ.get("SCDRS_GS", DATA / "psd_meta_magma_scdrs.gs")))
ONLY = os.environ.get("SCDRS_ATLAS")
print("gene sets:", {k: len(v[0]) for k, v in gs.items()}, flush=True)


def score_atlas(adata, atlas_key, group_cols):
    scdrs.preprocess(adata, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)
    rows = []
    for trait, (genes, weights) in gs.items():
        pres = [(g, w) for g, w in zip(genes, weights) if g in adata.var_names]
        if len(pres) < 50:
            print(f"  skip {trait}: {len(pres)} genes", flush=True); continue
        df = scdrs.score_cell(data=adata, gene_list=[p[0] for p in pres], gene_weight=[p[1] for p in pres], ctrl_match_key="mean_var",
                              n_ctrl=1000, weight_opt="vs", return_ctrl_raw_score=False, return_ctrl_norm_score=True, verbose=False).fillna(0.0)
        for gcol in group_cols:
            res = scdrs.method.downstream_group_analysis(adata=adata, df_full_score=df, group_cols=[gcol], fdr_thresholds=[0.05, 0.1, 0.2])[gcol]
            res = res.reset_index().rename(columns={"index": "group_level", gcol: "group_level", "group": "group_level"})
            res["atlas"], res["trait"], res["group_col"], res["n_genes_matched"] = atlas_key, trait, gcol, len(pres)
            rows.append(res)
        print(f"  {atlas_key} {trait}: {len(pres)} genes scored", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(RES / f"{TAG}_group_stats_{atlas_key}.tsv", sep="\t", index=False)
    return out


def load_gse287304():
    d = DATA / "GSE287304"; d.mkdir(exist_ok=True)
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287304/suppl/"
    for f in ("GSE287304_1_matrix.mtx.gz", "GSE287304_1_features.tsv.gz", "GSE287304_1_barcodes.tsv.gz"):
        if not (d / f).exists():
            print("  downloading", f, flush=True); urllib.request.urlretrieve(base + f, d / f)
    feats = pd.read_csv(d / "GSE287304_1_features.tsv.gz", sep="\t", header=None, names=["id", "name", "ftype"])
    bcs = pd.read_csv(d / "GSE287304_1_barcodes.tsv.gz", sep="\t", header=None, names=["barcode"])
    X = scipy.io.mmread(gzip.open(d / "GSE287304_1_matrix.mtx.gz")).tocsc().T.tocsr()      # cells x features
    rna = (feats.ftype == "Gene Expression").values; hto = (feats.ftype == "Antibody Capture").values
    Xr = X[:, rna]; keep = np.asarray((Xr > 0).sum(axis=1)).ravel() >= 500                  # same pre-filter as the original GSE287304 analysis
    Xr, Xh, bc = Xr[keep], X[keep][:, hto].toarray(), bcs.barcode.values[keep]
    names = feats.loc[hto, "id"].values
    h3 = Xh[:, [int(np.where(names == n)[0][0]) for n in ("HTO_C0252", "HTO_C0253", "HTO_C0254")]].astype(float)
    clr = np.log1p(h3); clr -= clr.mean(axis=0, keepdims=True)
    srt = np.sort(clr, axis=1); umi = np.asarray(Xr.sum(axis=1)).ravel()
    lab = np.array(["DSC (Dermal Sheath Cup)", "DP (Dermal Papilla)", "UDS (Upper Dermal Sheath)"])[np.argmax(clr, axis=1)]
    m = (srt[:, -1] > 0.5) & ((srt[:, -1] - srt[:, -2]) >= 0.35) & (umi >= 1000)
    a = ad.AnnData(X=Xr[m].astype(np.float32), obs=pd.DataFrame({"compartment": lab[m]}, index=bc[m]),
                   var=pd.DataFrame(index=feats.loc[rna, "name"].astype(str).values))
    a.var_names_make_unique(); sc.pp.filter_cells(a, min_genes=500)
    print("  GSE287304 singlets:", a.n_obs, a.obs.compartment.value_counts().to_dict(), flush=True)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a); sc.pp.pca(a, n_comps=30); sc.pp.neighbors(a)
    return a


for key, (url, desc) in s03.ATLAS_URLS.items():
    if ONLY and key != ONLY:
        continue
    h5 = DATA / f"{key}.h5ad"
    if not h5.exists():
        print(f"[download] {desc}", flush=True); urllib.request.urlretrieve(url, h5)
    print(f"[atlas] {desc}", flush=True)
    a = s03.prepare_adata_with_symbols(h5)
    score_atlas(a, key, s03.choose_group_columns(a)); del a
if ONLY in (None, "GSE274880_HS"):
    print("[atlas] GSE274880 HS + healthy skin", flush=True)
    a = s09.load_and_annotate_gse274880()
    score_atlas(a, "GSE274880_HS", ["cell_type", "cell_type_condition", "condition_group"]); del a
if ONLY in (None, "GSE287304_hair_follicle_mesenchyme"):
    print("[atlas] GSE287304 hair-follicle mesenchyme", flush=True)
    a = load_gse287304()
    score_atlas(a, "GSE287304_hair_follicle_mesenchyme", ["compartment"])
if ONLY is None:
    merge_atlases()
print("SCDRS_MAGMA_DONE")
