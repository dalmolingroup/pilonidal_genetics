#!/usr/bin/env python3
"""Step 06b: summary tables derived from the outputs of steps 01-05 (previously produced ad hoc).

  (1) LDSC h2 / intercepts and genetic correlations parsed from results/ldsc_logs (step 02b), with liability-scale h2
      -> results/ldsc_v2_heritability.tsv, results/ldsc_v2_genetic_correlations.tsv
  (2) lead p-values corrected for the LDSC intercept of each meta-analysis (chi2 / intercept)
      -> results/psd_meta_loci_intercept_corrected.tsv
  (3) inverse-variance meta-analysis of the non-FinnGen cohorts at each lead -> results/psd_meta_loci_nonfinngen_replication.tsv
  (4) LAVA bivariate results without MHC blocks, with Bonferroni and BH-FDR -> results/lava_bivariate_psd_noMHC.tsv
  (5) MHC flag on the trait colocalizations (results/coloc_psd_vs_traits.tsv, rewritten in place)
"""
import glob, os, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

BASE = Path(__file__).resolve().parents[1]
RES = BASE / "results"
LOGS = Path(os.environ.get("LDSC_LOG_DIR", RES / "ldsc_logs"))
# MHC: chr6 28,477,797-33,448,354 (GRCh38 GRC region), padded to 28-34 Mb
is_mhc = lambda chrom, pos: (str(chrom) == "6") and (28_000_000 <= int(pos) <= 34_000_000)

# ---------- (1) LDSC tables ----------
rows = []
for f in sorted(glob.glob(f"{LOGS}/h2_*.log")):
    t = open(f).read(); name = os.path.basename(f)[3:-4]
    g = lambda pat: re.search(pat, t)
    h2 = g(r"Total Observed scale h2: ([-\d.e]+) \(([-\d.e]+)\)"); lam = g(r"Lambda GC: ([\d.]+)"); mc = g(r"Mean Chi\^2: ([\d.]+)")
    ic = g(r"Intercept: ([\d.]+) \(([\d.]+)\)"); ra = g(r"Ratio: ([-\d.]+) \(([\d.]+)\)|Ratio < 0"); nsnp = g(r"After merging with regression SNP LD, (\d+) SNPs remain")
    rows.append({"trait": name, "n_snps": nsnp and int(nsnp.group(1)), "mean_chi2": mc and float(mc.group(1)), "lambda_gc": lam and float(lam.group(1)),
                 "h2_obs": h2 and float(h2.group(1)), "h2_obs_se": h2 and float(h2.group(2)), "h2_z": h2 and float(h2.group(1)) / float(h2.group(2)),
                 "intercept": ic and float(ic.group(1)), "intercept_se": ic and float(ic.group(2)),
                 "ratio": (ra.group(1) if ra and ra.group(1) else ("<0" if ra else None)), "ratio_se": ra and ra.group(2)})
h = pd.DataFrame(rows)

def liab(h2o, P, K):
    z = norm.pdf(norm.ppf(1 - K)); return h2o * K**2 * (1 - K)**2 / (P * (1 - P) * z**2)
P_samp = {"PSD_META_EUR": 6771 / (6771 + 1374212), "PSD_META_ALL": 7868 / (7868 + 1598492), "PSD_MVP_EUR": 1433 / 450253, "PSD_MVP_MultiAncestry": 2073 / 631604,
          "PSD_UKB_EUR": 1466 / 458440, "PSD_FinnGen_R13": 3872 / 472290}
for K in (0.005, 0.01):
    h[f"h2_liab_K{K}"] = [liab(r.h2_obs, P_samp[r.trait], K) if r.trait in P_samp else None for r in h.itertuples()]
    h[f"h2_liab_K{K}_se"] = [liab(r.h2_obs_se, P_samp[r.trait], K) if r.trait in P_samp else None for r in h.itertuples()]
h.to_csv(RES / "ldsc_v2_heritability.tsv", sep="\t", index=False, float_format="%.5g")

rg_rows = []
for f in sorted(glob.glob(f"{LOGS}/rg_*.log")):
    m = re.search(r"Summary of Genetic Correlation Results\n(.*?)\n\n", open(f).read(), re.S)
    if not m:
        print("no summary in", f); continue
    lines = [l.split() for l in m.group(1).strip().split("\n")]
    for l in lines[1:]:
        d = dict(zip(lines[0], l))
        d["p1"] = os.path.basename(d["p1"]).replace(".sumstats.gz", ""); d["p2"] = os.path.basename(d["p2"]).replace(".sumstats.gz", "")
        rg_rows.append(d)
rg = pd.DataFrame(rg_rows)
for c in ["rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]:
    if c in rg: rg[c] = pd.to_numeric(rg[c], errors="coerce")
rg.to_csv(RES / "ldsc_v2_genetic_correlations.tsv", sep="\t", index=False)
print(f"LDSC: {len(h)} h2 rows, {len(rg)} rg rows from {LOGS}", flush=True)

# ---------- (2) intercept-corrected lead p-values ----------
cl = pd.read_csv(RES / "psd_meta_locus_classification.tsv", sep="\t")
icpt = {m: float(h.loc[h.trait == f"PSD_META_{m}", "intercept"].iloc[0]) for m in ("ALL", "EUR")}
print("LDSC intercepts:", icpt, flush=True)
out = []
for r in cl.itertuples():
    best = None
    for m in ("ALL", "EUR"):
        p = getattr(r, f"p_META_{m}")
        if pd.notna(p):
            pc = chi2.sf(norm.isf(p / 2) ** 2 / icpt[m], 1)
            if best is None or pc < best[1]: best = (m, pc, p)
    out.append({"status": r.status, "chr": r.chr, "pos": r.pos, "gene": r.nearest_gene, "rsid": r.rsid, "best_meta": best[0], "p_raw": best[2],
                "p_intercept_corrected": best[1], "gws_after_correction": best[1] < 5e-8, "is_mhc": is_mhc(r.chr, r.pos)})
o = pd.DataFrame(out).sort_values("p_intercept_corrected")
o.to_csv(RES / "psd_meta_loci_intercept_corrected.tsv", sep="\t", index=False, float_format="%.3g")
x = o[~o.is_mhc]
print(f"GWS after intercept correction (non-MHC): {int(x.gws_after_correction.sum())}/{len(x)}", flush=True)

# ---------- (3) non-FinnGen check ----------
L = {m: pd.read_csv(RES / f"psd_meta_{m}_loci.tsv", sep="\t") for m in ("ALL", "EUR")}
rows = []
for r in cl.itertuples():
    src = None
    for m, cohorts in (("ALL", ["MVP_MultiAncestry", "UKB_EUR", "GH_SAS", "FinnGen_R13"]), ("EUR", ["MVP_EUR", "UKB_EUR", "FinnGen_R13"])):
        d = L[m]; hit = d[(d.chr.astype(str) == str(r.chr)) & ((d.pos - r.pos).abs() < 500_000) & (d.p_meta < 5e-8)]
        if len(hit): src = (m, cohorts, hit.sort_values("p_meta").iloc[0]); break
    m, cohorts, hit = src
    nonfg = [c for c in cohorts if c != "FinnGen_R13"]
    b = np.array([hit[f"{c}_beta"] for c in nonfg], float); se = np.array([hit[f"{c}_se"] for c in nonfg], float)
    ok = np.isfinite(b) & np.isfinite(se) & (se > 0)
    w = 1 / se[ok]**2; bm = (w * b[ok]).sum() / w.sum(); z = bm * np.sqrt(w.sum())
    same = np.sign(bm) == np.sign(hit.beta_meta)
    rows.append({"status": r.status, "locus": f"chr{r.chr}:{r.pos}", "gene": r.nearest_gene, "rsid": r.rsid, "meta": m, "p_meta": hit.p_meta,
                 "FinnGen_p": hit["FinnGen_R13_p"], "nonFinnGen_cohorts": "+".join(c for c, k in zip(nonfg, ok) if k), "nonFG_OR": np.exp(bm),
                 "nonFG_p": 2 * norm.sf(abs(z)), "nonFG_same_dir": bool(same), "nonFG_rep_p05": bool(same and 2 * norm.sf(abs(z)) < 0.05)})
t = pd.DataFrame(rows)
for c in ("p_meta", "FinnGen_p", "nonFG_p"): t[c] = t[c].map(lambda v: f"{v:.1e}" if pd.notna(v) else "")
t["nonFG_OR"] = t.nonFG_OR.round(3)
t.to_csv(RES / "psd_meta_loci_nonfinngen_replication.tsv", sep="\t", index=False)
nov = t[t.status == "NOVEL"]
print(f"NOVEL loci replicating outside FinnGen (same direction, p<0.05): {int(nov.nonFG_rep_p05.sum())}/{len(nov)}", flush=True)

# ---------- (4) LAVA without MHC ----------
def bh(p):
    p = np.asarray(p); n = len(p); o_ = np.argsort(p); r_ = np.empty(n); r_[o_] = np.arange(1, n + 1)
    q = p * n / r_; qs = np.minimum.accumulate(q[o_][::-1])[::-1]; res = np.empty(n); res[o_] = np.minimum(qs, 1); return res
lv = pd.read_csv(RES / "lava_bivariate_psd.tsv", sep="\t")
lv = lv[~((lv.chr == 6) & (lv.stop > 28_400_000) & (lv.start < 33_500_000))].copy()   # GRCh37 MHC blocks
lv["p_bonf"] = (lv.p * len(lv)).clip(upper=1); lv["fdr"] = bh(lv.p.values)
lv.to_csv(RES / "lava_bivariate_psd_noMHC.tsv", sep="\t", index=False)
print(f"LAVA non-MHC: {len(lv)} tests, FDR<0.05: {int((lv.fdr < 0.05).sum())}, Bonferroni<0.05: {int((lv.p_bonf < 0.05).sum())}", flush=True)

# ---------- (5) MHC flag on trait colocalizations ----------
co = pd.read_csv(RES / "coloc_psd_vs_traits.tsv", sep="\t")
co["is_mhc"] = [is_mhc(c, p) for c, p in zip(co.chr, co.pos)]
co.to_csv(RES / "coloc_psd_vs_traits.tsv", sep="\t", index=False)
print("DERIVED_TABLES_DONE", flush=True)
