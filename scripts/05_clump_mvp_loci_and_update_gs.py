#!/usr/bin/env python3
"""
Step 5: Distance-based locus clumping (500 kb window) of MVP Multi-Ancestry GWAS
(GCST90480439, 2,073 cases across EUR, AFR, AMR) and South Asian GWAS (GCST90727025),
stratified by common/low-frequency (EAF >= 0.005) vs ultra-rare variants, mapping lead
SNPs to nearest GRCh38 protein-coding genes via Ensembl REST API, and adding empirical
GWAS gene sets to data/psd_scdrs_genesets.gs.
"""
import json
import math
import urllib.request
from pathlib import Path
import pandas as pd
from scipy.stats import norm

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"


def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.loads(r.read().decode("utf-8"))


def map_region_to_genes(chrom: str, pos: int, window: int = 100000):
    start = max(1, pos - window)
    end = pos + window
    url = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}?feature=gene;content-type=application/json"
    try:
        genes = fetch_json(url)
        coding = [g for g in genes if g.get("biotype") == "protein_coding" and g.get("external_name")]
        if not coding:
            coding = [g for g in genes if g.get("external_name")]
        if not coding:
            return "", ""
        # sort by distance from pos to gene interval
        def dist(g):
            gs, ge = int(g["start"]), int(g["end"])
            if gs <= pos <= ge:
                return 0
            return min(abs(pos - gs), abs(pos - ge))
        coding.sort(key=dist)
        nearest = coding[0].get("external_name", "")
        all_names = ",".join(dict.fromkeys(g.get("external_name", "") for g in coding[:4]))
        return nearest, all_names
    except Exception:
        return "", ""


def clump_variants(df: pd.DataFrame, window_bp: int = 500000):
    df = df.sort_values("p_value", ascending=True).reset_index(drop=True)
    clumped = []
    for _, row in df.iterrows():
        c = str(row["chromosome"])
        p = int(row["base_pair_location"])
        if any(str(x["chromosome"]) == c and abs(int(x["base_pair_location"]) - p) < window_bp for x in clumped):
            continue
        clumped.append(row)
    return pd.DataFrame(clumped)


def main():
    mvp_path = RESULTS_DIR / "gwas_MVP_multiancestry_GCST90480439_top_variants.tsv"
    print(f"[1/3] Loading {mvp_path}...")
    df = pd.read_csv(mvp_path, sep="\t", dtype=str)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["base_pair_location"] = pd.to_numeric(df["base_pair_location"], errors="coerce")
    df["eaf"] = pd.to_numeric(df["effect_allele_frequency"], errors="coerce")
    df = df.dropna(subset=["p_value", "base_pair_location"]).copy()
    df["maf"] = df["eaf"].apply(lambda x: min(x, 1.0 - x) if pd.notnull(x) else 0.0)

    # 1. Common well-powered loci (MAF >= 2%, p < 5e-5)
    df_common = df[(df["maf"] >= 0.02) & (df["p_value"] < 5e-5)].copy()
    clumped_common = clump_variants(df_common, window_bp=500000)
    print(f"  Found {len(clumped_common)} independent common loci (MAF >= 2%, p < 5e-5) in MVP multi-ancestry GWAS")

    records = []
    mvp_gene_weights = {}
    for _, row in clumped_common.head(50).iterrows():
        chrom = str(row["chromosome"])
        pos = int(row["base_pair_location"])
        pval = float(row["p_value"])
        z = -norm.ppf(pval / 2.0)
        nearest, nearby = map_region_to_genes(chrom, pos, window=150000)
        if nearest:
            mvp_gene_weights[nearest] = max(mvp_gene_weights.get(nearest, 0.0), round(z, 4))
        records.append({
            "chrom": chrom,
            "pos_grch38": pos,
            "rsid": row.get("rsid", ""),
            "effect_allele": row.get("effect_allele", ""),
            "other_allele": row.get("other_allele", ""),
            "odds_ratio": row.get("odds_ratio", ""),
            "eaf": row.get("effect_allele_frequency", ""),
            "maf": round(float(row["maf"]), 5),
            "p_value": pval,
            "z_score": round(z, 3),
            "direction_EUR_AFR_AMR": str(row.get("direction", "")).strip(),
            "nearest_coding_gene": nearest,
            "genes_within_150kb": nearby,
        })

    df_loci = pd.DataFrame(records)
    out_loci = RESULTS_DIR / "mvp_multiancestry_independent_loci.tsv"
    df_loci.to_csv(out_loci, sep="\t", index=False)
    print(f"  Saved {len(df_loci)} annotated MVP independent loci -> {out_loci}")
    print(df_loci.head(20).to_string(index=False))

    # Append MVP_MultiAncestry_GWAS gene set to psd_scdrs_genesets.gs
    gs_path = DATA_DIR / "psd_scdrs_genesets.gs"
    if mvp_gene_weights and gs_path.exists():
        lines = [l for l in gs_path.read_text().splitlines() if not l.startswith("PSD_MVP_MultiAncestry_GWAS\t")]
        sorted_mvp = sorted(mvp_gene_weights.items(), key=lambda x: -x[1])
        lines.append("PSD_MVP_MultiAncestry_GWAS\t" + ",".join(f"{g}:{w:.4f}" for g, w in sorted_mvp))
        gs_path.write_text("\n".join(lines) + "\n")
        print(f"  Added PSD_MVP_MultiAncestry_GWAS ({len(sorted_mvp)} genes) to {gs_path}")


if __name__ == "__main__":
    main()
