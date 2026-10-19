#!/usr/bin/env python3
"""
Step 11d (lightweight, runs anywhere): merge META_ALL and META_EUR genome-wide significant loci into one locus list,
classify each as KNOWN (within 500 kb of a prior PSD association; data/psd_known_loci_prior_art.tsv) or NOVEL,
and record whether any single cohort reaches genome-wide significance on its own.
Inputs : results/psd_meta_{ALL,EUR}_loci.tsv, data/psd_known_loci_prior_art.tsv
Output : results/psd_meta_locus_classification.tsv
"""
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
WIN = 500_000
known = pd.read_csv(BASE / "data" / "psd_known_loci_prior_art.tsv", sep="\t")
known["chr"] = known["chr"].astype(str)

rows = []
for m in ("ALL", "EUR"):
    d = pd.read_csv(RES / f"psd_meta_{m}_loci.tsv", sep="\t")
    d = d[d.p_meta < 5e-8].copy(); d["meta"] = m; d["chr"] = d["chr"].astype(str)
    rows.append(d)
d = pd.concat(rows, ignore_index=True).sort_values("p_meta")
loci = []
for r in d.itertuples():
    hit = next((l for l in loci if l["chr"] == r.chr and abs(l["pos"] - int(r.pos)) < WIN), None)
    if hit is None:
        hit = {"chr": r.chr, "pos": int(r.pos), "rsid": r.rsid, "a1": r.a1, "a2": r.a2, "eaf": r.eaf_meta, "or": r.or_meta,
               "nearest_gene": r.nearest_coding_gene, "genes_250kb": r.genes_within_250kb, "p_META_ALL": None, "p_META_EUR": None}
        loci.append(hit)
    hit[f"p_META_{r.meta}"] = r.p_meta if hit[f"p_META_{r.meta}"] is None else min(hit[f"p_META_{r.meta}"], r.p_meta)
    for c in ("MVP_MultiAncestry", "MVP_EUR", "UKB_EUR", "FinnGen_R13", "GH_SAS"):
        if f"{c}_p" in d.columns and pd.notna(getattr(r, f"{c}_p", None)):
            hit[f"{c}_p"] = min(hit.get(f"{c}_p", 1.0), getattr(r, f"{c}_p"))
        if f"p_loo_{c}" in d.columns and pd.notna(getattr(r, f"p_loo_{c}", None)):
            hit[f"p_without_{c}"] = min(hit.get(f"p_without_{c}", 1.0), getattr(r, f"p_loo_{c}"))
    hit["i2_max"] = max(hit.get("i2_max", 0), r.i2 if pd.notna(r.i2) else 0)
    hit["direction_" + r.meta] = r.direction
out = pd.DataFrame(loci)
def classify(r):
    k = known[(known.chr == r.chr) & ((known.pos_b38 - r.pos).abs() < WIN)]
    return ("KNOWN", "; ".join(f"{a} [{b}]" for a, b in zip(k.locus_label, k.source))) if len(k) else ("NOVEL", "")
out[["status", "prior_art"]] = out.apply(lambda r: pd.Series(classify(r)), axis=1)
single = [c for c in out.columns if c.endswith("_p") and not c.startswith("p_")]
out["any_single_cohort_gws"] = (out[single] < 5e-8).any(axis=1)
out["is_mhc"] = (out.chr == "6") & out.pos.between(25_000_000, 34_000_000)
out = out.sort_values(["status", "p_META_ALL"], na_position="last")
out.to_csv(RES / "psd_meta_locus_classification.tsv", sep="\t", index=False, float_format="%.3g")
print(f"{len(out)} distinct GWS loci: {int((out.status=='KNOWN').sum())} known, {int((out.status=='NOVEL').sum())} novel; "
      f"novel loci not GWS in any single cohort: {int(((out.status=='NOVEL') & ~out.any_single_cohort_gws).sum())}")
pd.set_option("display.width", 250)
print(out[["status", "chr", "pos", "rsid", "nearest_gene", "eaf", "or", "p_META_ALL", "p_META_EUR"] + single + ["i2_max", "prior_art"]].to_string(index=False))
