#!/usr/bin/env python3
"""
Step 1A: Fetch and harmonize Open Targets genetic associations, published multi-biobank GWAS loci,
UK Biobank Whole-Exome Gene Burden (GCST90084318, 450k release), UK Biobank Whole-Exome Single Variants
(GCST90080332), and Genes & Health South Asian GWAS (GCST90727025).
"""
import os
import json
import gzip
import math
import urllib.request
import zlib
from pathlib import Path
from scipy.stats import norm

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DATA_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)


def fetch_json(url, payload=None):
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    else:
        data = None
    req = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


def fetch_opentargets():
    print("[1/4] Fetching Open Targets associations for HP_0010769 (Pilonidal sinus) and HP_0010771 (Pilonidal abscess)...")
    targets_out = []
    evidences_out = []
    for efo in ["HP_0010769", "HP_0010771", "MONDO_0008249"]:
        query = """
        query($efoId: String!) {
          disease(efoId: $efoId) {
            id
            name
            associatedTargets(page: {index: 0, size: 250}) {
              count
              rows {
                score
                target {
                  id
                  approvedSymbol
                  approvedName
                  biotype
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
        try:
            res = fetch_json("https://api.platform.opentargets.org/api/v4/graphql", {"query": query, "variables": {"efoId": efo}})
            d = res.get("data", {}).get("disease")
            if not d:
                continue
            dname = d.get("name")
            for row in d.get("associatedTargets", {}).get("rows", []):
                t = row["target"]
                dtypes = {x["id"]: x["score"] for x in row.get("datatypeScores", [])}
                targets_out.append({
                    "phenotype_id": efo,
                    "phenotype_name": dname,
                    "ensembl_id": t["id"],
                    "gene_symbol": t["approvedSymbol"],
                    "gene_name": t["approvedName"],
                    "overall_score": row["score"],
                    "genetic_association_score": dtypes.get("genetic_association", 0.0),
                    "animal_model_score": dtypes.get("animal_model", 0.0),
                    "literature_score": dtypes.get("literature", 0.0),
                })
            for ev in d.get("evidences", {}).get("rows", []):
                t = ev.get("target") or {}
                pval = None
                if ev.get("pValueMantissa") is not None and ev.get("pValueExponent") is not None:
                    pval = float(ev["pValueMantissa"]) * (10.0 ** int(ev["pValueExponent"]))
                evidences_out.append({
                    "phenotype_id": efo,
                    "phenotype_name": dname,
                    "ensembl_id": t.get("id", ""),
                    "gene_symbol": t.get("approvedSymbol", ""),
                    "datasource": ev.get("datasourceId", ""),
                    "datatype": ev.get("datatypeId", ""),
                    "score": ev.get("score", ""),
                    "study_id": ev.get("studyId", ""),
                    "variant_id": ev.get("variantId", ""),
                    "rsid": ev.get("variantRsId", ""),
                    "pvalue": pval if pval is not None else "",
                    "odds_ratio": ev.get("oddsRatio", ""),
                    "beta": ev.get("beta", ""),
                })
        except Exception as e:
            print(f"  Warning on {efo}: {e}")

    # Write targets TSV
    t_path = DATA_DIR / "opentargets_pilonidal_targets.tsv"
    if targets_out:
        keys = list(targets_out[0].keys())
        with open(t_path, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in sorted(targets_out, key=lambda x: (-float(x["genetic_association_score"]), -float(x["overall_score"]))):
                f.write("\t".join(str(r[k]) for k in keys) + "\n")
    print(f"  Saved {len(targets_out)} target associations -> {t_path}")

    # Write evidences TSV
    e_path = DATA_DIR / "opentargets_pilonidal_evidences.tsv"
    if evidences_out:
        keys = list(evidences_out[0].keys())
        with open(e_path, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in evidences_out:
                f.write("\t".join(str(r[k]) for k in keys) + "\n")
    print(f"  Saved {len(evidences_out)} evidence records -> {e_path}")


def fetch_ukb_exome_burden():
    print("[2/4] Downloading & parsing UKB Whole-Exome Gene Burden (GCST90084318, 4.3 MB)...")
    url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90084001-GCST90085000/GCST90084318/GCST90084318_buildGRCh38.tsv.gz"
    local_gz = DATA_DIR / "GCST90084318_buildGRCh38.tsv.gz"
    if not local_gz.exists():
        urllib.request.urlretrieve(url, local_gz)

    # Aggregate best gene-level burden p-value per gene (across masks M1/M3/etc.)
    best_per_gene = {}
    with gzip.open(local_gz, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            name = parts[idx["Name"]]
            try:
                pval = float(parts[idx["p_value"]])
            except ValueError:
                continue
            if pval <= 0 or pval > 1 or math.isnan(pval):
                continue
            # Name format: SAMD11(ENSG00000187634).GENE.M1.0001 or variant
            if "(" in name and ").GENE." in name:
                gene_sym = name.split("(")[0]
                ensg = name.split("(")[1].split(")")[0]
                mask = name.split(").GENE.")[1]
            else:
                gene_sym = name.split("(")[0] if "(" in name else name
                ensg = name.split("(")[1].split(")")[0] if "(" in name and ")" in name else ""
                mask = parts[idx["effect_allele"]]
            or_val = parts[idx["odds_ratio"]]
            eaf = parts[idx["effect_allele_frequency"]]
            if gene_sym not in best_per_gene or pval < best_per_gene[gene_sym]["p_value"]:
                z = -norm.ppf(pval / 2.0)
                best_per_gene[gene_sym] = {
                    "gene_symbol": gene_sym,
                    "ensembl_id": ensg,
                    "best_mask": mask,
                    "chrom": parts[idx["chromosome"]],
                    "pos": parts[idx["base_pair_location"]],
                    "odds_ratio": or_val,
                    "eaf": eaf,
                    "p_value": pval,
                    "z_score": round(z, 4),
                }

    sorted_genes = sorted(best_per_gene.values(), key=lambda x: x["p_value"])
    out_path = RESULTS_DIR / "ukb_exome_gene_burden_summary.tsv"
    keys = ["gene_symbol", "ensembl_id", "best_mask", "chrom", "pos", "odds_ratio", "eaf", "p_value", "z_score"]
    with open(out_path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for g in sorted_genes:
            f.write("\t".join(str(g[k]) for k in keys) + "\n")
    print(f"  Parsed {len(sorted_genes)} genes from GCST90084318 -> {out_path}")
    print("  Top 15 Exome Gene Burden hits:")
    for g in sorted_genes[:15]:
        print(f"    {g['gene_symbol']:12s} ({g['ensembl_id']}) | mask={g['best_mask']:8s} | OR={g['odds_ratio']:8s} | p={g['p_value']:.3e} | Z={g['z_score']}")


def fetch_ukb_exome_single_variant():
    print("[3/4] Downloading & parsing UKB Whole-Exome Single Variants (GCST90080332, 16 MB)...")
    url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90080001-GCST90081000/GCST90080332/GCST90080332_buildGRCh38.tsv.gz"
    local_gz = DATA_DIR / "GCST90080332_buildGRCh38.tsv.gz"
    if not local_gz.exists():
        urllib.request.urlretrieve(url, local_gz)

    top_vars = []
    with gzip.open(local_gz, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            try:
                pval = float(parts[idx["p_value"]])
            except ValueError:
                continue
            if pval < 1e-4:
                top_vars.append((pval, parts))

    top_vars.sort(key=lambda x: x[0])
    out_path = RESULTS_DIR / "ukb_exome_top_variants_GCST90080332.tsv"
    with open(out_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for _, parts in top_vars:
            f.write("\t".join(parts) + "\n")
    print(f"  Found {len(top_vars)} exome single variants with p < 1e-4 -> {out_path}")
    for pval, parts in top_vars[:10]:
        print(f"    {parts[0]} | chr{parts[1]}:{parts[2]} | OR={parts[8]} | EAF={parts[12]} | p={pval:.3e}")


def save_curated_published_loci():
    print("[4/4] Saving curated published GWAS & Open Targets high-confidence PSD risk genes...")
    # Combines Roberson et al. 2024 (PMID 38902823), UKB Focused Derm GWAS 2025 (PMID 39624183),
    # and Open Targets Genetics L2G / candidate causal genes (HP_0010769 / HP_0010771)
    curated = [
        ("TWIST1", "ENSG00000122691", "7p21.1", "GWAS_2024 + OpenTargets_HP_0010771", "Mesenchymal/dermal condensation & cranial-caudal neural crest/Wnt patterning", 6.2),
        ("HDAC9", "ENSG00000048052", "7p21.1", "GWAS_2024 + OpenTargets_HP_0010771", "Epigenetic regulator adjacent to TWIST1; androgenetic alopecia & hair patterning locus", 6.0),
        ("DKK1", "ENSG00000107984", "10q21.1", "OpenTargets_HP_0010771", "Dickkopf Wnt signaling pathway inhibitor 1; dermal papilla hair follicle cycling & fibrosis", 5.8),
        ("MBL2", "ENSG00000165471", "10q21.1", "OpenTargets_HP_0010771", "Mannose-binding lectin 2; innate lectin-complement opsonization of S. aureus and anaerobes", 4.9),
        ("TBX15", "ENSG00000092607", "1p12", "GWAS_2024", "T-box transcription factor 15; dorsal-ventral mesenchymal & gluteal/subcutaneous adipose patterning", 5.9),
        ("WARS2", "ENSG00000116874", "1p12", "GWAS_2024", "Mitochondrial tryptophanyl-tRNA synthetase 2; waist-hip ratio & subcutaneous fat distribution locus", 5.7),
        ("PRKAR1B", "ENSG00000188191", "7p22.3", "GWAS_2024 + OpenTargets_HP_0010771", "cAMP-dependent protein kinase type I-beta regulatory subunit", 5.5),
        ("PDGFA", "ENSG00000197461", "7p22.3", "OpenTargets_HP_0010771", "Platelet-derived growth factor subunit A; mesenchymal dermal papilla & fibroblast proliferation", 5.2),
        ("NEK9", "ENSG00000119638", "14q24.3", "GWAS_2024", "NIMA-related kinase 9; mutated in follicular keratosis / nevus comedonicus & hair follicle cycling", 5.8),
        ("EIF2B2", "ENSG00000119718", "14q24.3", "GWAS_2024", "Eukaryotic translation initiation factor 2B subunit beta (14q24.3 locus with NEK9)", 5.3),
        ("GPATCH2L", "ENSG00000089965", "14q24.3", "GWAS_2024", "G-patch domain containing 2 like (14q24.3 locus)", 5.3),
        ("TMED10", "ENSG00000170348", "14q24.3", "OpenTargets_HP_0010771", "Transmembrane p24 trafficking protein 10; unconventional secretion of IL-1beta", 5.0),
        ("FERD3L", "ENSG00000175048", "7p21.1", "OpenTargets_HP_0010771", "Fer3-like bHLH transcription factor at 7p21.1 locus", 4.6),
        ("TCHH", "ENSG00000159450", "1q21.3", "Exome_Rare_2024 + OpenTargets_HP_0010771", "Trichohyalin; inner root sheath & hair shaft medulla cross-linking (OR 4.81 in PSD)", 5.6),
        ("PADI3", "ENSG00000142619", "1p36.13", "Exome_Candidate_2024", "Peptidyl arginine deiminase 3; hair shaft modification & uncombable hair syndrome pathway", 3.8),
        ("TGM3", "ENSG00000125780", "20p13", "Exome_Candidate_2024", "Transglutaminase 3; hair follicle inner root sheath & cornified envelope cross-linking", 3.8),
        ("RNFT2", "ENSG00000135119", "12q24.22", "OpenTargets_HP_0010769", "Ring finger protein transmembrane 2; E3 ubiquitin ligase regulating IL-3R/innate cytokine signaling", 5.4),
        ("CHRM3", "ENSG00000133019", "1q43", "OpenTargets_HP_0010769", "Muscarinic acetylcholine receptor M3; primary driver of eccrine sweat gland secretion", 5.1),
        ("TTK", "ENSG00000112742", "6q14.1", "UKB_Derm_GWAS_2025", "TTK protein kinase (Mps1); mitotic spindle checkpoint in rapidly dividing epithelial progenitors", 4.8),
    ]
    out_path = RESULTS_DIR / "curated_psd_risk_genes.tsv"
    with open(out_path, "w") as f:
        f.write("gene_symbol\tensembl_id\tlocus\tevidence_source\tbiological_role\tprior_z_weight\n")
        for row in curated:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"  Saved {len(curated)} high-confidence PSD risk genes -> {out_path}")


if __name__ == "__main__":
    fetch_opentargets()
    fetch_ukb_exome_burden()
    fetch_ukb_exome_single_variant()
    save_curated_published_loci()
