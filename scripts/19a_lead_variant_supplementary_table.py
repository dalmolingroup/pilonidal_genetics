#!/usr/bin/env python3
"""Step 19a: supplementary table of the 25 non-MHC PSD lead variants, with full statistical reporting in both meta-analyses.

For each lead variant (results/psd_locus_to_gene.tsv; GRCh38) and each meta-analysis (META_EUR, META_ALL):
  - GRCh38 and GRCh37 positions (GRCh37 from the Ensembl GRCh37 cache written by step 17a), rsID, effect/other allele, EAF
  - fixed-effect IVW OR (95% CI) and p, as in step 11b
  - DerSimonian-Laird random-effects OR (95% CI) and p, tau^2, Cochran's Q, heterogeneity p and I^2
  - per-cohort OR (95% CI) and p
Indels keep the VCF position (the padding base before the variant), so they sit 1 bp before the Ensembl start.
Output: results/psd_lead_variants_supplementary.tsv
"""
import json, os
from pathlib import Path
import duckdb, numpy as np, pandas as pd
from scipy.stats import norm, chi2

BASE = Path(__file__).resolve().parents[1]
COHORTS = {"EUR": ["MVP_EUR", "UKB_EUR", "FinnGen_R13"], "ALL": ["MVP_MultiAncestry", "UKB_EUR", "FinnGen_R13", "GH_SAS"]}

l2g = pd.read_csv(BASE / "results/psd_locus_to_gene.tsv", sep="\t")
loc = l2g.locus.str.extract(r"chr(\w+):(\d+)").astype(int)
leads = pd.DataFrame({"locus": l2g.locus, "chr": loc[0], "pos": loc[1], "lead_rsid": l2g.lead_rsid, "nearest_gene": l2g.nearest_gene,
                      "status": l2g.status})
cls = pd.read_csv(BASE / "results/psd_meta_locus_classification.tsv", sep="\t")[["chr", "pos", "a1", "a2"]]
leads = leads.merge(cls, on=["chr", "pos"], how="left")
b37 = json.loads((BASE / "data/processed/psd_leads_grch37.json").read_text())
con = duckdb.connect()
con.execute(f"SET threads={os.environ.get('NCPU', 4)}")
con.register("leads", leads[["chr", "pos"]])

def ci(b, se):
    return f"{np.exp(b):.3f} ({np.exp(b - 1.96 * se):.3f}-{np.exp(b + 1.96 * se):.3f})"

# per-cohort effects at the lead positions, with the same MAF filter and within-cohort de-duplication as step 11b
coh = []
for c in sorted({c for v in COHORTS.values() for c in v}):
    f = BASE / f"data/processed/harmonized/PSD_{c}.harmonized.tsv.gz"
    coh.append(con.execute(f"""
        WITH raw AS (
          SELECT TRY_CAST(chr AS INTEGER) AS chr, CAST(pos AS BIGINT) AS pos, a1, a2, CAST(beta AS DOUBLE) AS beta,
                 CAST(se AS DOUBLE) AS se, TRY_CAST(p AS DOUBLE) AS p
          FROM read_csv('{f}', delim='\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false,
                        ignore_errors=true, max_line_size=1000000)
          WHERE CAST(eaf AS DOUBLE) BETWEEN 0.01 AND 0.99 AND CAST(se AS DOUBLE) > 0)
        SELECT '{c}' AS cohort, r.chr, r.pos, r.a1, r.a2, arg_min(r.beta, r.se) AS beta, MIN(r.se) AS se, arg_min(r.p, r.se) AS p
        FROM raw r JOIN leads l ON r.chr = l.chr AND r.pos = l.pos
        GROUP BY r.chr, r.pos, r.a1, r.a2""").df())
coh = pd.concat(coh)

rows = []
for meta, cohorts in COHORTS.items():
    f = BASE / f"data/processed/meta/psd_meta_{meta}.tsv.gz"
    d = con.execute(f"SELECT m.* FROM read_csv_auto('{f}', delim='\t', header=true) m JOIN leads USING (chr, pos)").df()
    for c in cohorts:                                                                     # attach cohort columns, matched on alleles
        x = coh[coh.cohort == c][["chr", "pos", "a1", "a2", "beta", "se", "p"]].rename(columns={"beta": f"{c}_beta", "se": f"{c}_se", "p": f"{c}_p"})
        d = d.merge(x, on=["chr", "pos", "a1", "a2"], how="left")
    for r in leads.itertuples():
        hit = d[(d.chr == r.chr) & (d.pos == r.pos)]
        if len(hit) > 1:                                                                  # multi-allelic site: keep the reported lead allele
            hit = hit[(hit.a1 == r.a1) & (hit.a2 == r.a2)]
        out = {"locus": r.locus, "status": r.status, "nearest_gene": r.nearest_gene, "rsid": r.lead_rsid, "chr": r.chr,
               "pos_grch38": r.pos, "pos_grch37": b37.get(r.lead_rsid), "meta": f"META_{meta}"}
        if hit.empty:
            rows.append({**out, "note": "variant not in this meta-analysis"}); continue
        h = hit.iloc[0]
        b = np.array([h[f"{c}_beta"] for c in cohorts], float); se = np.array([h[f"{c}_se"] for c in cohorts], float)
        ok = ~np.isnan(b) & ~np.isnan(se); b, se = b[ok], se[ok]; k = len(b)
        w = 1 / se**2; bf = (w * b).sum() / w.sum(); sef = w.sum() ** -0.5
        Q = (w * (b - bf) ** 2).sum(); tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum())) if k > 1 else 0.0
        wr = 1 / (se**2 + tau2); br = (wr * b).sum() / wr.sum(); ser = wr.sum() ** -0.5
        out.update(effect_allele=h.a1, other_allele=h.a2, eaf=round(h.eaf_meta, 4), k_cohorts=k, n_eff=round(h.n_eff),
                   or_fixed_95ci=ci(bf, sef), p_fixed=2 * norm.sf(abs(bf / sef)),
                   or_random_95ci=ci(br, ser), p_random=2 * norm.sf(abs(br / ser)), tau2=tau2,
                   cochran_q=Q, p_het=chi2.sf(Q, k - 1) if k > 1 else np.nan, i2=max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0,
                   check_beta_vs_11b=round(bf - h.beta_meta, 6))
        for c in cohorts:
            bc, sc = h[f"{c}_beta"], h[f"{c}_se"]
            out[f"{c}_or_95ci"] = "" if pd.isna(bc) else ci(bc, sc)
            out[f"{c}_p"] = h[f"{c}_p"]
        rows.append(out)

res = pd.DataFrame(rows)
res.to_csv(BASE / "results/psd_lead_variants_supplementary.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(res[["locus", "meta", "rsid", "pos_grch37", "or_fixed_95ci", "p_fixed", "or_random_95ci", "p_random", "i2", "p_het"]].to_string(index=False))
print("max |beta_fixed - beta_11b|:", res.check_beta_vs_11b.abs().max())
print("LEAD_TABLE_DONE")
