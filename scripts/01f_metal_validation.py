#!/usr/bin/env python3
"""Step 01f: reproduce both meta-analyses of step 01b with METAL (2011-03-25) and compare them variant by variant.

METAL is given the per-cohort records that 01b combines (same QC: MAF >= 0.01, exact duplicates collapsed to the
smallest SE, multi-allelic sites excluded within and across cohorts), with markers named chr:pos so that METAL does
its own allele matching. Analysis: SCHEME STDERR (fixed-effect inverse-variance weighting), GENOMICCONTROL OFF,
ANALYZE HETEROGENEITY.
Inputs : data/processed/harmonized/PSD_<cohort>.harmonized.tsv.gz, data/processed/meta/psd_meta_{ALL,EUR}.tsv.gz,
         results/psd_meta_{ALL,EUR}_loci.tsv (01b)
Outputs: data/processed/metal/META_{ALL,EUR}1.tbl, results/psd_meta_metal_concordance.tsv,
         results/psd_meta_metal_leads.tsv
"""
import os, subprocess, time
from pathlib import Path
import duckdb
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
HARM, META = BASE / "data/processed/harmonized", BASE / "data/processed/meta"
OUT, RES = BASE / "data/processed/metal", BASE / "results"
OUT.mkdir(parents=True, exist_ok=True)
METAL = os.environ.get("METAL", "metal")
METAS = {"ALL": ["MVP_MultiAncestry", "UKB_EUR", "FinnGen_R13", "GH_SAS"], "EUR": ["MVP_EUR", "UKB_EUR", "FinnGen_R13"]}
MAF_MIN, GWS = 0.01, 5e-8

con = duckdb.connect()
con.execute(f"SET memory_limit='{os.environ.get('DUCKDB_MEM', '12GB')}'; SET threads={os.environ.get('DUCKDB_THREADS', '8')};")
con.execute(f"SET temp_directory='{os.environ.get('DUCKDB_TMP', OUT / 'duckdb_tmp')}';")

# ---- per-cohort records after the 01b QC (same SQL as 01b, steps 1 and the multi-allelic filter) ----
parts = [f"""SELECT '{c}' AS cohort, chr, CAST(pos AS BIGINT) AS pos, a1, a2, CAST(beta AS DOUBLE) AS beta,
                    CAST(se AS DOUBLE) AS se, CAST(eaf AS DOUBLE) AS eaf
             FROM read_csv('{HARM}/PSD_{c}.harmonized.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='',
                           comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
             WHERE CAST(se AS DOUBLE) > 0 AND isfinite(CAST(beta AS DOUBLE))""" for c in sorted({c for v in METAS.values() for c in v})]
con.execute("CREATE TABLE long_dedup AS SELECT cohort, chr, pos, a1, a2, arg_min(beta, se) AS beta, MIN(se) AS se, arg_min(eaf, se) AS eaf "
            "FROM (" + " UNION ALL ".join(parts) + ") GROUP BY cohort, chr, pos, a1, a2")
con.execute(f"""CREATE TABLE long AS SELECT * FROM long_dedup
                ANTI JOIN (SELECT cohort, chr, pos FROM long_dedup GROUP BY cohort, chr, pos HAVING COUNT(*) > 1) USING (cohort, chr, pos)
                WHERE eaf BETWEEN {MAF_MIN} AND {1 - MAF_MIN}""")
con.execute("DROP TABLE long_dedup")
print(time.strftime("%H:%M:%S"), "per-cohort QC done", flush=True)

summary, leads_out = [], []
for m, cohorts in METAS.items():
    # 01b output for this meta-analysis; its (chr, pos) set already excludes positions with different allele pairs across cohorts
    con.execute(f"""CREATE OR REPLACE TABLE ours AS SELECT CAST(chr AS VARCHAR) AS chr, pos, a1, a2, k, beta_meta, se_meta, p_meta, i2
                    FROM read_csv('{META}/psd_meta_{m}.tsv.gz', delim='\\t', header=true, types={{'chr': 'VARCHAR'}})""")
    files = []
    for c in cohorts:
        f = OUT / f"{m}_{c}.tsv"
        con.execute(f"""COPY (SELECT chr || ':' || pos AS MARKER, a1 AS A1, a2 AS A2, beta AS BETA, se AS SE, eaf AS EAF
                              FROM long WHERE cohort = '{c}' AND (chr, pos) IN (SELECT chr, pos FROM ours))
                        TO '{f}' (DELIMITER '\\t', HEADER)""")
        files.append(f)
    script = OUT / f"META_{m}.metal"
    script.write_text("\n".join(["SCHEME STDERR", "GENOMICCONTROL OFF", "AVERAGEFREQ ON", "SEPARATOR TAB", "MARKER MARKER",
                                 "ALLELE A1 A2", "EFFECT BETA", "STDERR SE", "FREQ EAF"] + [f"PROCESS {f}" for f in files]
                                + [f"OUTFILE {OUT}/META_{m} .tbl", "ANALYZE HETEROGENEITY", "QUIT", ""]))
    print(time.strftime("%H:%M:%S"), f"running METAL for META_{m}", flush=True)
    with open(OUT / f"META_{m}.log", "w") as log:
        subprocess.run([METAL, str(script)], stdout=log, stderr=subprocess.STDOUT, check=True)
    for f in files:
        f.unlink()

    # METAL reports the effect for Allele1 (lower case); align it to our a1
    con.execute(f"""CREATE OR REPLACE TABLE metal AS
        SELECT split_part(MarkerName, ':', 1) AS chr, CAST(split_part(MarkerName, ':', 2) AS BIGINT) AS pos,
               upper(Allele1) AS al1, upper(Allele2) AS al2, Effect, StdErr, "P-value" AS p, HetISq,
               length(replace(Direction, '?', '')) AS k
        FROM read_csv('{OUT}/META_{m}1.tbl', delim='\\t', header=true)""")
    con.execute("""CREATE OR REPLACE TABLE cmp AS
        SELECT o.*, CASE WHEN x.al1 = o.a1 THEN x.Effect ELSE -x.Effect END AS beta_metal, x.StdErr AS se_metal,
               x.p AS p_metal, x.HetISq AS i2_metal, x.k AS k_metal
        FROM ours o JOIN metal x ON o.chr = x.chr AND o.pos = x.pos
             AND ((x.al1 = o.a1 AND x.al2 = o.a2) OR (x.al1 = o.a2 AND x.al2 = o.a1))""")
    s = con.execute(f"""SELECT COUNT(*) AS n_ours, (SELECT COUNT(*) FROM metal) AS n_metal, (SELECT COUNT(*) FROM cmp) AS n_matched,
               (SELECT COUNT(*) FROM cmp WHERE k <> k_metal) AS n_k_differ,
               (SELECT corr(beta_meta, beta_metal) FROM cmp) AS r_beta,
               (SELECT corr(-log10(p_meta), -log10(p_metal)) FROM cmp WHERE p_meta > 0 AND p_metal > 0) AS r_log10p,
               (SELECT MAX(abs(beta_meta - beta_metal)) FROM cmp) AS max_abs_d_beta,
               (SELECT MAX(abs(se_meta - se_metal)) FROM cmp) AS max_abs_d_se,
               (SELECT MAX(abs(log10(p_meta) - log10(p_metal))) FROM cmp WHERE p_meta > 0 AND p_metal > 0) AS max_abs_d_log10p,
               (SELECT MAX(abs(i2 - i2_metal)) FROM cmp WHERE k >= 2) AS max_abs_d_i2,
               (SELECT COUNT(*) FROM cmp WHERE k >= 2 AND p_meta < {GWS}) AS n_gws_ours,
               (SELECT COUNT(*) FROM cmp WHERE k_metal >= 2 AND p_metal < {GWS}) AS n_gws_metal,
               (SELECT COUNT(*) FROM cmp WHERE k >= 2 AND (p_meta < {GWS}) <> (p_metal < {GWS})) AS n_gws_discordant
        FROM ours""").df().assign(meta=f"META_{m}")
    summary.append(s)
    print(s.T.to_string(), flush=True)

    lead = pd.read_csv(RES / f"psd_meta_{m}_loci.tsv", sep="\t", dtype={"chr": str})
    lead = lead[lead.p_meta < GWS][["chr", "pos", "a1", "a2", "rsid"]]
    con.register("lead", lead)
    leads_out.append(con.execute(f"""SELECT 'META_{m}' AS meta, l.rsid, c.chr, c.pos, c.a1, c.a2, c.beta_meta, c.beta_metal, c.se_meta, c.se_metal,
                                            c.p_meta, c.p_metal, c.i2, c.i2_metal
                                     FROM lead l JOIN cmp c ON c.chr = l.chr AND c.pos = l.pos AND c.a1 = l.a1 AND c.a2 = l.a2""").df())

pd.concat(summary)[["meta"] + list(summary[0].columns[:-1])].to_csv(RES / "psd_meta_metal_concordance.tsv", sep="\t", index=False, float_format="%.6g")
pd.concat(leads_out).to_csv(RES / "psd_meta_metal_leads.tsv", sep="\t", index=False, float_format="%.6g")
print("METAL_VALIDATION_DONE", flush=True)
