#!/usr/bin/env python3
"""
Step 2: Stream/parse GWAS Catalog summary statistics (Genes & Health SAS GCST90727025
and MVP Multi-Ancestry GCST90480439), map top variants to GRCh38 genes, and generate
multi-threshold & pathway-stratified .gs files in official scDRS format.
"""
import os
import gzip
import math
import json
import urllib.request
from pathlib import Path
from scipy.stats import norm

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"


def stream_gcst90727025():
    """Download/scan Genes & Health South Asian GWAS (GCST90727025, 148 MB) for p < 1e-4."""
    url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90727001-GCST90728000/GCST90727025/GCST90727025.tsv"
    out_top = RESULTS_DIR / "gwas_SAS_GCST90727025_top_variants.tsv"
    if out_top.exists():
        print(f"[1/3] SAS top variants already cached at {out_top}")
        return

    print("[1/3] Streaming Genes & Health South Asian GWAS (GCST90727025, 457 cases / 43,569 controls)...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    top_rows = []
    with urllib.request.urlopen(req, timeout=60) as r:
        header = r.readline().decode("utf-8").strip().split("\t")
        idx = {k: i for i, k in enumerate(header)}
        p_idx = idx["p_value"]
        for raw_line in r:
            line = raw_line.decode("utf-8", errors="ignore").rstrip("\n")
            parts = line.split("\t")
            if len(parts) <= p_idx:
                continue
            try:
                pval = float(parts[p_idx])
            except ValueError:
                continue
            if pval < 1e-4:
                top_rows.append((pval, parts))

    top_rows.sort(key=lambda x: x[0])
    with open(out_top, "w") as f:
        f.write("\t".join(header) + "\n")
        for _, parts in top_rows:
            f.write("\t".join(parts) + "\n")
    print(f"  Saved {len(top_rows)} SAS GWAS variants (p < 1e-4) -> {out_top}")
    for pval, parts in top_rows[:10]:
        print(f"    chr{parts[0]}:{parts[1]} ({parts[8]}) | beta={parts[4]} | EAF={parts[6]} | p={pval:.3e}")


def stream_gcst90480439():
    """Stream MVP Multi-Ancestry GWAS (GCST90480439, 2,073 cases EUR+AFR+AMR, 834 MB gz) for p < 1e-4."""
    url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480439/GCST90480439.tsv.gz"
    out_top = RESULTS_DIR / "gwas_MVP_multiancestry_GCST90480439_top_variants.tsv"
    if out_top.exists():
        print(f"[2/3] MVP multi-ancestry top variants already cached at {out_top}")
        return

    print("[2/3] Streaming MVP Multi-Ancestry GWAS (GCST90480439, 2,073 cases / 629,531 controls)...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    top_rows = []
    with urllib.request.urlopen(req, timeout=120) as r:
        with gzip.GzipFile(fileobj=r) as gz:
            header = gz.readline().decode("utf-8").strip().split("\t")
            idx = {k: i for i, k in enumerate(header)}
            p_idx = idx["p_value"]
            for raw_line in gz:
                line = raw_line.decode("utf-8", errors="ignore").rstrip("\n")
                parts = line.split("\t")
                if len(parts) <= p_idx:
                    continue
                p_str = parts[p_idx]
                if p_str == "#NA" or "e-" not in p_str and not p_str.startswith("0.000"):
                    # fast pre-filter before float conversion
                    try:
                        pval = float(p_str)
                    except ValueError:
                        continue
                    if pval >= 1e-4:
                        continue
                else:
                    try:
                        pval = float(p_str)
                    except ValueError:
                        continue
                    if pval >= 1e-4:
                        continue
                top_rows.append((pval, parts))

    top_rows.sort(key=lambda x: x[0])
    with open(out_top, "w") as f:
        f.write("\t".join(header) + "\n")
        for _, parts in top_rows:
            f.write("\t".join(parts) + "\n")
    print(f"  Saved {len(top_rows)} MVP multi-ancestry variants (p < 1e-4) -> {out_top}")
    for pval, parts in top_rows[:15]:
        print(f"    chr{parts[0]}:{parts[1]} ({parts[8]}) | OR={parts[4]} | dir={parts[20]} | p={pval:.3e}")


def build_scdrs_genesets():
    """Build official scDRS .gs files combining Locus-Prioritized, OpenTargets, and Exome Gene Burden sets."""
    print("[3/3] Building scDRS gene set (.gs) files...")
    # 1. Load curated high-confidence GWAS + OpenTargets + Rare coding loci
    curated_genes = {}
    with open(RESULTS_DIR / "curated_psd_risk_genes.tsv") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            curated_genes[parts[0]] = float(parts[5])

    # 2. Load Open Targets genetic association targets (HP_0010769 & HP_0010771)
    ot_genes = {}
    ot_path = DATA_DIR / "opentargets_pilonidal_targets.tsv"
    if ot_path.exists():
        with open(ot_path) as f:
            header = f.readline().strip().split("\t")
            idx = {k: i for i, k in enumerate(header)}
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < len(header):
                    continue
                sym = parts[idx["gene_symbol"]]
                gscore = float(parts[idx["genetic_association_score"]])
                oscore = float(parts[idx["overall_score"]])
                if gscore > 0:
                    z_equiv = round(2.5 + 10.0 * gscore, 4)
                    ot_genes[sym] = max(ot_genes.get(sym, 0.0), z_equiv)

    # 3. Load UKB Whole-Exome Gene Burden (GCST90084318)
    exome_genes = []
    with open(RESULTS_DIR / "ukb_exome_gene_burden_summary.tsv") as f:
        header = f.readline().strip().split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            sym = parts[idx["gene_symbol"]]
            z = float(parts[idx["z_score"]])
            pval = float(parts[idx["p_value"]])
            exome_genes.append((sym, z, pval))

    # Construct scDRS traits dictionary: trait_name -> dict of {gene: weight}
    traits = {}

    # Set A: High-confidence GWAS + OpenTargets L2G + Rare coding (19 core genes)
    set_a = dict(curated_genes)
    for g, w in ot_genes.items():
        if g not in set_a:
            set_a[g] = w
    traits["PSD_GWAS_Locus_Prioritized"] = set_a

    # Set B: UKB Exome Rare Coding Burden Top 50 & Top 200
    traits["PSD_Exome_Burden_Top50"] = {sym: z for sym, z, _ in exome_genes[:50]}
    traits["PSD_Exome_Burden_Top200"] = {sym: z for sym, z, _ in exome_genes[:200]}

    # Set C: Integrated Polygenic + Exome + GWAS (Top 100 combined)
    combined = dict(set_a)
    for sym, z, _ in exome_genes[:100]:
        combined[sym] = max(combined.get(sym, 0.0), z)
    traits["PSD_Integrated_Top100"] = combined

    gs_path = DATA_DIR / "psd_scdrs_genesets.gs"
    with open(gs_path, "w") as f:
        f.write("TRAIT\tGENESET\n")
        for tname, gdict in traits.items():
            sorted_items = sorted(gdict.items(), key=lambda x: -x[1])
            gstr = ",".join(f"{g}:{w:.4f}" for g, w in sorted_items)
            f.write(f"{tname}\t{gstr}\n")
    print(f"  Wrote {len(traits)} scDRS gene sets -> {gs_path}")


if __name__ == "__main__":
    stream_gcst90727025()
    build_scdrs_genesets()
    stream_gcst90480439()
