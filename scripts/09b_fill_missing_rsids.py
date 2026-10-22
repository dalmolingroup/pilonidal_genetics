#!/usr/bin/env python3
"""Step 09b: fill rsIDs missing from the result tables (indels that lacked an rsID in the harmonised cohort files).

Missing variants are resolved with the Ensembl VEP REST endpoint (GRCh38), which left-aligns/normalises the VCF-style
indel and returns the co-located dbSNP rsID. Results are cached in data/processed/indel_rsids_grch38.json.
Every missing variant in these tables is an indel; one position (chr20:22130748, +A and +AA) is a multi-allelic repeat
and both alleles share one rsID. Run after any rerun of steps 01b-01e, 04b or 05a, which regenerate these tables.
"""
import json, urllib.request
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
R = BASE / "results"
CACHE = BASE / "data/processed/indel_rsids_grch38.json"
# table -> list of (rsid column, chr column, position column)
TABLES = {
    "psd_meta_EUR_loci": [("rsid", "chr", "pos")], "psd_meta_ALL_loci": [("rsid", "chr", "pos")],
    "psd_meta_locus_classification": [("rsid", "chr", "pos")], "psd_meta_loci_intercept_corrected": [("rsid", "chr", "pos")],
    "psd_meta_loci_nonfinngen_replication": [("rsid", None, "locus")], "psd_finemap_credible_sets": [("rsid", "chr", "pos")],
    "psd_finemap_vep_annotation": [("rsid", "chr", "pos")], "psd_finemap_vep_by_rsid": [("rsid", "chr", "pos")],
    "psd_finemap_locus_summary": [("abf_top_rsid", "chr", "abf_top_pos")],
    "coloc_psd_vs_traits": [("psd_lead", "chr", "pos")], "eqtl_coloc_psd_gtex": [("psd_lead", "chr", "pos")],
    "psd_published_loci_replication": [("META_ALL_best_rsid", "chr", "META_ALL_best_pos"), ("META_EUR_best_rsid", "chr", "META_EUR_best_pos")],
}

tabs = {t: pd.read_csv(R / f"{t}.tsv", sep="\t", dtype=str, keep_default_na=False, na_values=[""]) for t in TABLES}   # text kept as-is ("#NA" etc.)
# variants to resolve: alleles come from the tables that carry them
need = {}
for t in ["psd_meta_EUR_loci", "psd_meta_ALL_loci", "psd_finemap_credible_sets", "psd_finemap_vep_annotation"]:
    d = tabs[t]
    for r in d[d.rsid.isna()].itertuples():
        need[f"{r.chr} {r.pos} . {r.a1} {r.a2} . . ."] = (r.chr, int(r.pos))
cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
todo = [v for v in need if v not in cache]
if todo:
    req = urllib.request.Request("https://rest.ensembl.org/vep/human/region", data=json.dumps({"variants": todo}).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    for x in json.load(urllib.request.urlopen(req, timeout=120)):
        key = next(v for v in todo if v.startswith(x["input"].split(" .")[0] + " .") and v == x["input"])
        ids = [c["id"] for c in x.get("colocated_variants", []) if c["id"].startswith("rs")]
        cache[key] = ids[0] if len(ids) == 1 else None
    CACHE.write_text(json.dumps(cache, indent=1))
by_pos = {}
for v, (c, p) in need.items():
    if cache.get(v):
        assert by_pos.get((c, p), cache[v]) == cache[v], f"two rsIDs at {c}:{p}"
        by_pos[(c, p)] = cache[v]
unresolved = [v for v in need if not cache.get(v)]

# the HDAC9 locus summary row has no lead position column; its lead is the abf_top variant (same position)
for t, specs in TABLES.items():
    d = tabs[t]; n = 0
    for col, cchr, cpos in specs:
        miss = d[col].isna() & d[cpos].notna()
        if cchr is None:                                                                  # position only as "chrC:POS"
            cp = d.loc[miss, cpos].str.extract(r"chr(\w+):(\d+)"); pairs = zip(cp[0], cp[1])
        else:
            pairs = zip(d.loc[miss, cchr], d.loc[miss, cpos])
        fill = [by_pos.get((c, int(float(p)))) for c, p in pairs]
        d.loc[miss, col] = fill; n += sum(f is not None for f in fill)
    if t == "psd_finemap_locus_summary":
        m = d.lead_rsid.isna() & d.locus.str.match(r"chr\w+_\d+_")
        pos = d.loc[m, "locus"].str.extract(r"chr(\w+)_(\d+)_")
        d.loc[m, "lead_rsid"] = [by_pos.get((c, int(p))) for c, p in zip(pos[0], pos[1])]; n += int(m.sum())
    d.to_csv(R / f"{t}.tsv", sep="\t", index=False)
    print(f"{t}: filled {n}")
print("resolved:", {f"{c}:{p}": r for (c, p), r in sorted(by_pos.items())})
print("unresolved:", unresolved)
