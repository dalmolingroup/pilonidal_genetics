#!/usr/bin/env python3
"""
Step 11b: Multi-cohort inverse-variance-weighted fixed-effects meta-analysis of Pilonidal Sinus Disease GWAS.

Cohorts (all GRCh38, harmonized by scripts/11a_harmonize_psd_cohorts.sh):
  MVP_MultiAncestry  GCST90480439  2,073 cases / 629,531 controls (EUR+AFR+AMR; MVP-internal meta)
  MVP_EUR            GCST90478766  1,433 cases / 448,820 controls (subset of MVP_MultiAncestry -> never combined with it)
  UKB_EUR            GCST90473908  1,466 cases / 456,974 controls (UKB WGS, non-Finnish European)
  FinnGen_R13        L12_PILONIDALCYST 3,872 cases / 468,418 controls
  GH_SAS             GCST90727025    457 cases /  43,569 controls (Genes & Health, British Bangladeshi/Pakistani)

Two independent-sample meta-analyses:
  META_ALL (multi-ancestry): MVP_MultiAncestry + UKB_EUR + FinnGen_R13 + GH_SAS   (7,868 cases)
  META_EUR (European only): MVP_EUR + UKB_EUR + FinnGen_R13                      (6,771 cases; used for EUR-LD fine-mapping)

Per-variant inclusion: MAF >= 0.01 within cohort. Discovery claims require >= 2 contributing cohorts.
Outputs: data/processed/meta/psd_meta_{ALL,EUR}.tsv.gz, results/psd_meta_*_loci.tsv, results/psd_meta_cohort_qc.tsv,
         results/psd_published_loci_replication.tsv, figures/psd_meta_manhattan_qq.png
"""
import json
import os
import math
import time
import urllib.request
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2 as chi2_dist

BASE = Path(__file__).resolve().parent.parent
HARM = BASE / "data" / "processed" / "harmonized"
META_DIR = BASE / "data" / "processed" / "meta"
RESULTS = BASE / "results"
FIGS = BASE / "figures"
for d in (META_DIR, RESULTS, FIGS):
    d.mkdir(parents=True, exist_ok=True)

COHORTS = {
    "MVP_MultiAncestry": dict(cases=2073, controls=629531, n_max=631604, ancestry="EUR+AFR+AMR"),
    "MVP_EUR":           dict(cases=1433, controls=448820, n_max=450253, ancestry="EUR"),
    "UKB_EUR":           dict(cases=1466, controls=456974, n_max=458440, ancestry="EUR"),
    "FinnGen_R13":       dict(cases=3872, controls=468418, n_max=472290, ancestry="FIN"),
    "GH_SAS":            dict(cases=457,  controls=43569,  n_max=44026,  ancestry="SAS"),
}
for c in COHORTS.values():
    c["n_eff"] = 4.0 / (1.0 / c["cases"] + 1.0 / c["controls"])

METAS = {
    "ALL": ["MVP_MultiAncestry", "UKB_EUR", "FinnGen_R13", "GH_SAS"],
    "EUR": ["MVP_EUR", "UKB_EUR", "FinnGen_R13"],
}
MAF_MIN = 0.01
GWS = 5e-8
SUGG = 1e-6
CLUMP_BP = 500_000

# Published PSD loci (Wasserman-type 2024 DCR GWAS: UKB+FinnGen+PennMed EUR, 2,835 cases; JID Innov 2025 UKB focused GWAS)
PUBLISHED_LOCI = {
    "HDAC9/TWIST1 (7p21.1)": ["HDAC9", "TWIST1"],
    "TBX15/WARS2 (1p12)": ["TBX15", "WARS2"],
    "PRKAR1B (7p22.3)": ["PRKAR1B"],
    "NEK9/GPATCH2L/EIF2B2 (14q24.3)": ["NEK9", "GPATCH2L", "EIF2B2"],
    "TTK (6q14.1; JID Innov 2025, rare coding)": ["TTK"],
}


def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def ensembl_gene_coords(symbol):
    try:
        g = fetch_json(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json")
        return str(g["seq_region_name"]), int(g["start"]), int(g["end"])
    except Exception:
        return None


def nearest_genes(chrom, pos, window=250_000):
    try:
        genes = fetch_json(f"https://rest.ensembl.org/overlap/region/human/{chrom}:{max(1,pos-window)}-{pos+window}?feature=gene;content-type=application/json")
    except Exception:
        return "", ""
    coding = [g for g in genes if g.get("biotype") == "protein_coding" and g.get("external_name")]
    def dist(g):
        s, e = int(g["start"]), int(g["end"])
        return 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
    coding.sort(key=dist)
    near = coding[0]["external_name"] if coding else ""
    within = ",".join(dict.fromkeys(g["external_name"] for g in coding[:6]))
    return near, within


def main():
    con = duckdb.connect()
    con.execute(f"SET memory_limit='{os.environ.get('DUCKDB_MEM','12GB')}'; SET threads={os.environ.get('DUCKDB_THREADS','12')};")
    con.execute(f"SET temp_directory='{META_DIR}/duckdb_tmp';")

    # ---------- 1. Load harmonized cohorts into one long table ----------
    print("[1/6] Loading harmonized cohorts...", flush=True)
    union_parts = []
    for name, meta in COHORTS.items():
        f = HARM / f"PSD_{name}.harmonized.tsv.gz"
        if not f.exists():
            raise SystemExit(f"missing harmonized file {f}")
        union_parts.append(f"""
            SELECT '{name}' AS cohort, chr, CAST(pos AS BIGINT) AS pos, a1, a2,
                   CAST(beta AS DOUBLE) AS beta, CAST(se AS DOUBLE) AS se, TRY_CAST(p AS DOUBLE) AS p,
                   CAST(eaf AS DOUBLE) AS eaf, TRY_CAST(n AS DOUBLE) AS n,
                   CASE WHEN rsid LIKE 'rs%' THEN split_part(rsid,'|',1) ELSE NULL END AS rsid,
                   CASE WHEN rsid LIKE '%|%' THEN split_part(rsid,'|',2) ELSE NULL END AS mvp_direction
            FROM read_csv('{f}', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
            WHERE CAST(eaf AS DOUBLE) BETWEEN {MAF_MIN} AND {1-MAF_MIN} AND CAST(se AS DOUBLE) > 0 AND isfinite(CAST(beta AS DOUBLE))
        """)
    con.execute("CREATE TABLE long_raw AS " + " UNION ALL ".join(union_parts))
    # De-duplicate within cohort: one row per (cohort, chr, pos, a1, a2), keeping the smallest-SE record
    # (multi-allelic splits / duplicated indels otherwise inflate k and double-count a cohort).
    con.execute("""
        CREATE TABLE long AS
        SELECT cohort, chr, pos, a1, a2,
               arg_min(beta, se) AS beta, MIN(se) AS se, arg_min(p, se) AS p, arg_min(eaf, se) AS eaf,
               arg_min(n, se) AS n, arg_min(rsid, se) AS rsid, arg_min(mvp_direction, se) AS mvp_direction
        FROM long_raw GROUP BY cohort, chr, pos, a1, a2
    """)
    n_raw, n_dedup = con.execute("SELECT (SELECT COUNT(*) FROM long_raw), (SELECT COUNT(*) FROM long)").fetchone()
    print(f"   de-duplicated within-cohort keys: {n_raw:,} -> {n_dedup:,} rows", flush=True)
    con.execute("DROP TABLE long_raw")

    # ---------- 2. Per-cohort QC ----------
    print("[2/6] Per-cohort QC (lambda_GC, mean chi2)...", flush=True)
    qc = con.execute("""
        SELECT cohort, COUNT(*) AS n_variants_maf01,
               AVG((beta/se)*(beta/se)) AS mean_chi2,
               quantile_cont((beta/se)*(beta/se), 0.5)/0.4549364 AS lambda_gc,
               SUM(CASE WHEN p < 5e-8 THEN 1 ELSE 0 END) AS n_gws_variants
        FROM long GROUP BY cohort ORDER BY cohort
    """).df()
    for name, meta in COHORTS.items():
        qc.loc[qc.cohort == name, ["n_cases", "n_controls", "n_eff", "ancestry"]] = [meta["cases"], meta["controls"], round(meta["n_eff"], 1), meta["ancestry"]]
    qc.to_csv(RESULTS / "psd_meta_cohort_qc.tsv", sep="\t", index=False)
    print(qc.to_string(index=False), flush=True)

    all_meta_lead = {}
    for meta_name, cohorts in METAS.items():
        print(f"[3/6] IVW fixed-effects meta-analysis: META_{meta_name} ({' + '.join(cohorts)})", flush=True)
        cohort_list = ",".join(f"'{c}'" for c in cohorts)
        order_case = " ".join(f"WHEN '{c}' THEN {i}" for i, c in enumerate(cohorts))
        dir_cols = ",\n".join(
            f"MAX(CASE WHEN cohort='{c}' THEN (CASE WHEN beta>0 THEN '+' ELSE '-' END) END) AS d{i}" for i, c in enumerate(cohorts)
        )
        neff_scale = " ".join(
            f"WHEN '{c}' THEN {COHORTS[c]['n_eff']} * LEAST(1.0, COALESCE(n,{COHORTS[c]['n_max']})/{COHORTS[c]['n_max']})" for c in cohorts
        )
        con.execute(f"DROP TABLE IF EXISTS meta_{meta_name}")
        con.execute(f"""
            CREATE TABLE meta_{meta_name} AS
            WITH l AS (
                SELECT *, 1.0/(se*se) AS w, (CASE cohort {neff_scale} END) AS n_eff_c
                FROM long WHERE cohort IN ({cohort_list})
            ),
            g AS (
                SELECT chr, pos, a1, a2,
                       COUNT(*) AS k,
                       SUM(w*beta)/SUM(w) AS beta_meta,
                       1.0/SQRT(SUM(w)) AS se_meta,
                       SUM(w*beta*beta) - (SUM(w*beta)*SUM(w*beta))/SUM(w) AS q_het,
                       SUM(n_eff_c*eaf)/SUM(n_eff_c) AS eaf_meta,
                       SUM(n_eff_c) AS n_eff,
                       MIN(p) AS min_cohort_p,
                       ANY_VALUE(rsid) AS rsid,
                       ANY_VALUE(mvp_direction) AS mvp_direction,
                       {dir_cols}
                FROM l GROUP BY chr, pos, a1, a2
            )
            SELECT chr, pos, a1, a2, rsid, k, beta_meta, se_meta, beta_meta/se_meta AS z_meta,
                   eaf_meta, n_eff, q_het, min_cohort_p, mvp_direction,
                   {"||".join(f"COALESCE(d{i},'?')" for i in range(len(cohorts)))} AS direction
            FROM g
        """)
        # p-values & heterogeneity in pandas (vectorized, accurate tails)
        df = con.execute(f"SELECT * FROM meta_{meta_name}").df()
        z = df["z_meta"].values
        df["p_meta"] = 2.0 * norm.sf(np.abs(z))
        logp = -(norm.logsf(np.abs(z)) + math.log(2.0)) / math.log(10.0)
        df["neglog10p_meta"] = logp
        dfree = (df["k"] - 1).clip(lower=0)
        df["p_het"] = np.where(dfree > 0, chi2_dist.sf(df["q_het"].clip(lower=0), dfree.clip(lower=1)), np.nan)
        df["i2"] = np.where(df["q_het"] > 0, np.clip((df["q_het"] - dfree) / df["q_het"], 0, 1) * 100, 0.0)
        df["or_meta"] = np.exp(df["beta_meta"])
        n_total = len(df)
        n_k2 = int((df["k"] >= 2).sum())
        print(f"   {n_total:,} variants; {n_k2:,} present in >=2 cohorts; genome-wide significant (k>=2): {int(((df.p_meta<GWS)&(df.k>=2)).sum()):,}; (any k): {int((df.p_meta<GWS).sum()):,}", flush=True)
        # meta-level inflation on k>=2
        sub = df[df.k >= 2]
        chi = sub["z_meta"].values ** 2
        print(f"   META_{meta_name} k>=2: mean chi2 = {chi.mean():.4f}, lambda_GC = {np.median(chi)/0.4549364:.4f}", flush=True)
        qc_row = pd.DataFrame([{"cohort": f"META_{meta_name}", "n_variants_maf01": len(sub), "mean_chi2": chi.mean(), "lambda_gc": np.median(chi)/0.4549364,
                                "n_gws_variants": int((sub.p_meta < GWS).sum()), "n_cases": sum(COHORTS[c]["cases"] for c in cohorts),
                                "n_controls": sum(COHORTS[c]["controls"] for c in cohorts), "n_eff": round(sum(COHORTS[c]["n_eff"] for c in cohorts), 1),
                                "ancestry": "+".join(sorted(set(COHORTS[c]["ancestry"] for c in cohorts)))}])
        qc_row.to_csv(RESULTS / "psd_meta_cohort_qc.tsv", sep="\t", index=False, mode="a", header=False)

        out_cols = ["chr", "pos", "a1", "a2", "rsid", "k", "direction", "mvp_direction", "beta_meta", "se_meta", "or_meta", "z_meta", "p_meta", "neglog10p_meta",
                    "eaf_meta", "n_eff", "q_het", "p_het", "i2", "min_cohort_p"]
        df = df[out_cols].sort_values(["chr", "pos"], key=lambda s: s.astype(int) if s.name == "chr" else s)
        out = META_DIR / f"psd_meta_{meta_name}.tsv.gz"
        df.to_csv(out, sep="\t", index=False, float_format="%.6g")
        print(f"   wrote {out}", flush=True)

        # ---------- 4. Lead loci (distance clumping) ----------
        print(f"[4/6] Clumping META_{meta_name} loci (p < {SUGG}, k >= 2, +/-{CLUMP_BP//1000} kb)...", flush=True)
        cand = df[(df.p_meta < SUGG) & (df.k >= 2)].sort_values("p_meta").copy()
        cand["chr_i"] = cand["chr"].astype(int)
        leads = []
        for _, r in cand.iterrows():
            if any(l["chr_i"] == r["chr_i"] and abs(int(l["pos"]) - int(r["pos"])) < CLUMP_BP for l in leads):
                continue
            leads.append(r)
        leads = pd.DataFrame(leads)
        if len(leads):
            # per-cohort stats for lead variants
            keys = ", ".join(f"('{r.chr}',{int(r.pos)},'{r.a1}','{r.a2}')" for r in leads.itertuples())
            per = con.execute(f"""
                SELECT cohort, chr, pos, a1, a2, beta, se, p, eaf FROM long
                WHERE (chr,pos,a1,a2) IN ({keys}) AND cohort IN ({cohort_list})
            """).df()
            per["pos"] = per["pos"].astype(int)
            rows = []
            for r in leads.itertuples():
                sub = per[(per.chr == r.chr) & (per.pos == int(r.pos)) & (per.a1 == r.a1) & (per.a2 == r.a2)]
                d = {c: v for c, v in r._asdict().items() if c not in ("Index", "chr_i")}
                for c in cohorts:
                    s = sub[sub.cohort == c]
                    if len(s):
                        s = s.iloc[0]
                        d[f"{c}_beta"] = round(float(s.beta), 4); d[f"{c}_se"] = round(float(s.se), 4); d[f"{c}_p"] = float(s.p) if pd.notna(s.p) else None; d[f"{c}_eaf"] = round(float(s.eaf), 4)
                    else:
                        d[f"{c}_beta"] = None; d[f"{c}_se"] = None; d[f"{c}_p"] = None; d[f"{c}_eaf"] = None
                # leave-one-cohort-out IVW p-values (robustness: is the locus supported without the largest cohort?)
                for c_out in cohorts:
                    bb, ww = [], []
                    for c in cohorts:
                        if c != c_out and d.get(f"{c}_se"):
                            bb.append(d[f"{c}_beta"]); ww.append(1.0 / d[f"{c}_se"] ** 2)
                    d[f"p_loo_{c_out}"] = float(2 * norm.sf(abs(np.sum(np.array(ww) * np.array(bb)) / np.sqrt(np.sum(ww))))) if ww else None
                d["n_cohorts_nominal_same_dir"] = int(sum(1 for c in cohorts if d.get(f"{c}_p") is not None and d[f"{c}_p"] < 0.05 and np.sign(d[f"{c}_beta"]) == np.sign(r.beta_meta)))
                d["significance"] = "genome-wide" if r.p_meta < GWS else "suggestive"
                near, within = nearest_genes(r.chr, int(r.pos))
                time.sleep(0.15)
                d["nearest_coding_gene"] = near; d["genes_within_250kb"] = within
                rows.append(d)
            leads = pd.DataFrame(rows)
        all_meta_lead[meta_name] = leads
        out_l = RESULTS / f"psd_meta_{meta_name}_loci.tsv"
        leads.to_csv(out_l, sep="\t", index=False, float_format="%.4g")
        print(f"   {len(leads)} independent loci (GWS: {int((leads.p_meta < GWS).sum()) if len(leads) else 0}) -> {out_l}", flush=True)
        if len(leads):
            show = leads[["chr", "pos", "rsid", "a1", "a2", "eaf_meta", "or_meta", "p_meta", "k", "direction", "i2", "n_cohorts_nominal_same_dir", "nearest_coding_gene", "significance"]]
            print(show.to_string(index=False), flush=True)

    # ---------- 5. Replication of published loci ----------
    print("[5/6] Replication look-up of published PSD loci in each cohort and meta...", flush=True)
    rep_rows = []
    for locus, genes in PUBLISHED_LOCI.items():
        coords = [ensembl_gene_coords(g) for g in genes]
        coords = [c for c in coords if c]
        if not coords:
            continue
        chrom = coords[0][0]; start = min(c[1] for c in coords) - CLUMP_BP; end = max(c[2] for c in coords) + CLUMP_BP
        row = {"published_locus": locus, "chr": chrom, "window_start": start, "window_end": end}
        for meta_name in METAS:
            r = con.execute(f"""SELECT pos, a1, a2, rsid, beta_meta, se_meta, k, direction, eaf_meta FROM meta_{meta_name}
                               WHERE chr='{chrom}' AND pos BETWEEN {start} AND {end} AND k >= 2
                               ORDER BY abs(beta_meta/se_meta) DESC LIMIT 1""").df()
            if len(r):
                r = r.iloc[0]; z = r.beta_meta / r.se_meta
                row[f"META_{meta_name}_best_rsid"] = r.rsid; row[f"META_{meta_name}_best_pos"] = int(r.pos); row[f"META_{meta_name}_best_p"] = 2 * norm.sf(abs(z))
                row[f"META_{meta_name}_best_or"] = math.exp(r.beta_meta); row[f"META_{meta_name}_best_dir"] = r.direction; row[f"META_{meta_name}_best_eaf"] = r.eaf_meta
        for c in COHORTS:
            r = con.execute(f"""SELECT pos, rsid, beta, se, p FROM long WHERE cohort='{c}' AND chr='{chrom}' AND pos BETWEEN {start} AND {end}
                               ORDER BY p ASC NULLS LAST LIMIT 1""").df()
            if len(r):
                r = r.iloc[0]; row[f"{c}_best_p"] = r.p; row[f"{c}_best_rsid"] = r.rsid; row[f"{c}_best_or"] = math.exp(r.beta)
                row[f"{c}_n_tests"] = int(con.execute(f"SELECT COUNT(*) FROM long WHERE cohort='{c}' AND chr='{chrom}' AND pos BETWEEN {start} AND {end}").fetchone()[0])
        rep_rows.append(row)
    rep = pd.DataFrame(rep_rows)
    rep.to_csv(RESULTS / "psd_published_loci_replication.tsv", sep="\t", index=False, float_format="%.4g")
    print(rep[[c for c in rep.columns if c.endswith("_best_p") or c == "published_locus"]].to_string(index=False), flush=True)

    # ---------- 6. Manhattan + QQ ----------
    print("[6/6] Plotting Manhattan & QQ...", flush=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(18, 10), gridspec_kw={"width_ratios": [3.2, 1]})
    for i, meta_name in enumerate(METAS):
        df = pd.read_csv(META_DIR / f"psd_meta_{meta_name}.tsv.gz", sep="\t", usecols=["chr", "pos", "p_meta", "neglog10p_meta", "k"])
        df = df[df.k >= 2]
        df["chr"] = df["chr"].astype(int)
        df = df.sort_values(["chr", "pos"])
        offsets = {}; off = 0; ticks = []
        for c in range(1, 23):
            sub = df[df.chr == c]
            offsets[c] = off
            ticks.append(off + (sub.pos.max() if len(sub) else 0) / 2)
            off += (sub.pos.max() if len(sub) else 0) + 5_000_000
        x = df.pos.values + df.chr.map(offsets).values
        ax = axes[i, 0]
        plot_mask = df.neglog10p_meta.values > 1.0
        ax.scatter(x[plot_mask], df.neglog10p_meta.values[plot_mask], s=2, c=np.where(df.chr.values[plot_mask] % 2 == 0, "#1f4e79", "#7fa6c9"), rasterized=True)
        ax.axhline(-math.log10(GWS), color="red", lw=0.8, ls="--"); ax.axhline(-math.log10(SUGG), color="orange", lw=0.6, ls=":")
        ax.set_xticks(ticks); ax.set_xticklabels(range(1, 23), fontsize=8); ax.set_ylabel("-log10(P)")
        cases = sum(COHORTS[c]["cases"] for c in METAS[meta_name])
        ax.set_title(f"PSD META_{meta_name}: {' + '.join(METAS[meta_name])} ({cases:,} cases), MAF>=1%, >=2 cohorts")
        leads = all_meta_lead[meta_name]
        if len(leads):
            for r in leads[leads.p_meta < GWS].itertuples():
                ax.annotate(r.nearest_coding_gene or r.rsid, (int(r.pos) + offsets[int(r.chr)], r.neglog10p_meta), fontsize=8, ha="center", va="bottom")
        # QQ
        p = np.sort(df.p_meta.values)
        exp = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        obs = -np.log10(np.clip(p, 1e-300, 1))
        axq = axes[i, 1]
        thin = np.concatenate([np.arange(0, min(len(p), 50000)), np.arange(50000, len(p), max(1, len(p) // 200000))]) if len(p) > 50000 else np.arange(len(p))
        axq.scatter(exp[::-1][thin], obs[::-1][thin], s=3, c="#1f4e79", rasterized=True)
        axq.plot([0, exp.max()], [0, exp.max()], c="red", lw=0.8)
        chi = norm.isf(p / 2) ** 2
        axq.set_title(f"QQ META_{meta_name} (lambda_GC={np.median(chi)/0.4549364:.3f})"); axq.set_xlabel("Expected -log10(P)"); axq.set_ylabel("Observed -log10(P)")
    plt.tight_layout()
    fig.savefig(FIGS / "psd_meta_manhattan_qq.png", dpi=160)
    print(f"   saved {FIGS / 'psd_meta_manhattan_qq.png'}", flush=True)
    con.close()


if __name__ == "__main__":
    main()
