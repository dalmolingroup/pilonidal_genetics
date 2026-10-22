#!/usr/bin/env python3
"""
Step 04a: full-density, rsID-keyed summary statistics for colocalization, LAVA and Genomic SEM
(cluster job; DuckDB). One file per trait in data/processed/shared/<TRAIT>.tsv.gz with columns
SNP CHR POS A1 A2 BETA SE Z P N EAF  (A1 = effect allele; MAF >= 1%; rsIDs present once only;
palindromic A/T, C/G SNPs with 0.4 < EAF < 0.6 removed). POS is GRCh38 except BMI/WHRadjBMI (GRCh37);
all downstream matching is by rsID + alleles, so builds are never mixed.

Traits:
  PSD_META_EUR / PSD_META_ALL : this study (scripts 01a-01b); N = total N scaled by per-SNP effective-N fraction
  HS        : GCST90726609 (MVP+UKB+FinnGen; 3,941 cases / 1,435,603 controls; 73% EUR)
  HS_FG     : FinnGen R13 L12_HIDRADENITISSUP (1,453 / 474,896; EUR-only sensitivity)
  ACNE_FG   : FinnGen R13 L12_ACNE (hospital-diagnosed; 4,812 / 474,896)
  ACNE_TL   : GCST90245818 Teder-Laving 2024 (34,422 / 364,991)
  BALDP1    : GCST90043616 UKB males, balding PATTERN 1 (= no balding; 66,172 / 140,864). Negative rg = positive with balding.
  BMI, WHRadjBMI : Pulit 2019 GIANT+UKB
Also writes results/hs_gws_loci.tsv (distance-clumped HS GCST90726609 loci, p < 5e-8, 500 kb).
"""
import os
from pathlib import Path
import duckdb

BASE = Path(__file__).resolve().parent.parent
S = BASE / "data" / "raw" / "sumstats"
META = BASE / "data" / "processed" / "meta"
OUT = BASE / "data" / "processed" / "shared"
OUT.mkdir(parents=True, exist_ok=True)
con = duckdb.connect()
con.execute(f"SET memory_limit='{os.environ.get('DUCKDB_MEM','48GB')}'; SET threads={os.environ.get('DUCKDB_THREADS','16')};")
con.execute(f"SET temp_directory='{OUT}/duckdb_tmp';")
O = "delim='\\t', header=true, all_varchar=true, quote='', escape='', strict_mode=false, ignore_errors=true, max_line_size=1000000"


def write(name, sql):
    out = OUT / f"{name}.tsv.gz"
    if out.exists() and out.stat().st_size > 1_000_000:
        print(f"{name}: exists, skipping", flush=True); return
    con.execute(f"""
      COPY (
        WITH raw AS ({sql}),
        clean AS (
          SELECT SNP, CHR, CAST(POS AS BIGINT) POS, upper(A1) A1, upper(A2) A2, BETA, SE, BETA/SE AS Z, P, CAST(N AS BIGINT) N, EAF FROM raw
          WHERE SNP LIKE 'rs%' AND SE > 0 AND isfinite(BETA) AND P IS NOT NULL AND N > 0 AND EAF BETWEEN 0.01 AND 0.99
            AND regexp_matches(upper(A1), '^[ACGT]+$') AND regexp_matches(upper(A2), '^[ACGT]+$') AND upper(A1) <> upper(A2)
            AND NOT ((upper(A1)||upper(A2)) IN ('AT','TA','CG','GC') AND EAF BETWEEN 0.4 AND 0.6)
        )
        SELECT * FROM clean WHERE SNP IN (SELECT SNP FROM clean GROUP BY SNP HAVING COUNT(*) = 1)
      ) TO '{out}' (DELIMITER '\\t', HEADER true, COMPRESSION gzip)""")
    n = con.execute(f"SELECT COUNT(*) FROM read_csv('{out}', delim='\\t', header=true)").fetchone()[0]
    print(f"{name}: {n:,} SNPs -> {out}", flush=True)


for m, ntot, neff in (("EUR", 1381000, 26920.0), ("ALL", 1606360, 31280.1)):
    write(f"PSD_META_{m}", f"""
      SELECT split_part(split_part(rsid, '|', 1), ',', 1) SNP, chr CHR, pos POS, a1 A1, a2 A2, CAST(beta_meta AS DOUBLE) BETA, CAST(se_meta AS DOUBLE) SE,
             CAST(p_meta AS DOUBLE) P, ROUND({ntot} * CAST(n_eff AS DOUBLE) / {neff}) N, CAST(eaf_meta AS DOUBLE) EAF
      FROM read_csv('{META}/psd_meta_{m}.tsv.gz', {O}) WHERE CAST(k AS INT) >= 2""")
write("HS", f"""
  SELECT rsid SNP, chromosome CHR, base_pair_location POS, effect_allele A1, other_allele A2, TRY_CAST(beta AS DOUBLE) BETA,
         TRY_CAST(standard_error AS DOUBLE) SE, TRY_CAST(p_value AS DOUBLE) P, 1439544 N, TRY_CAST(effect_allele_frequency AS DOUBLE) EAF
  FROM read_csv('{S}/HS_MVP_UKB_FinnGen_GCST90726609.h.tsv.gz', {O})""")
for name, f, n in (("HS_FG", "HS_FinnGen_R13_L12_HIDRADENITISSUP.gz", 1453 + 474896), ("ACNE_FG", "ACNE_FinnGen_R13_L12_ACNE.gz", 4812 + 474896)):
    write(name, f"""
      SELECT split_part(rsids, ',', 1) SNP, "#chrom" CHR, pos POS, alt A1, ref A2, TRY_CAST(beta AS DOUBLE) BETA, TRY_CAST(sebeta AS DOUBLE) SE,
             TRY_CAST(pval AS DOUBLE) P, {n} N, TRY_CAST(af_alt AS DOUBLE) EAF
      FROM read_csv('{S}/{f}', {O})""")
write("ACNE_TL", f"""
  SELECT rsid SNP, chromosome CHR, base_pair_location POS, effect_allele A1, other_allele A2, TRY_CAST(beta AS DOUBLE) BETA,
         TRY_CAST(standard_error AS DOUBLE) SE, TRY_CAST(p_value AS DOUBLE) P, 399413 N, TRY_CAST(effect_allele_frequency AS DOUBLE) EAF
  FROM read_csv('{S}/ACNE_TederLaving2024_GCST90245818.h.tsv.gz', {O})""")
write("BALDP1", f"""
  SELECT hm_rsid SNP, hm_chrom CHR, hm_pos POS, hm_effect_allele A1, hm_other_allele A2, TRY_CAST(hm_beta AS DOUBLE) BETA,
         TRY_CAST(standard_error AS DOUBLE) SE, TRY_CAST(p_value AS DOUBLE) P, COALESCE(TRY_CAST(n AS DOUBLE), 207036) N,
         TRY_CAST(hm_effect_allele_frequency AS DOUBLE) EAF
  FROM read_csv('{S}/AGA_UKB_GCST90043616.h.tsv.gz', {O})""")
for name, f in (("BMI", "BMI_GIANT_UKB_Pulit2019.txt.gz"), ("WHRadjBMI", "WHRadjBMI_GIANT_UKB_Pulit2019.txt.gz")):
    write(name, f"""
      SELECT split_part(SNP, ':', 1) SNP, CHR, POS, Tested_Allele A1, Other_Allele A2, TRY_CAST(BETA AS DOUBLE) BETA, TRY_CAST(SE AS DOUBLE) SE,
             TRY_CAST(P AS DOUBLE) P, TRY_CAST(N AS DOUBLE) N, TRY_CAST(Freq_Tested_Allele AS DOUBLE) EAF
      FROM read_csv('{S}/{f}', {O.replace("delim='\\t'", "delim=' '")}) WHERE COALESCE(TRY_CAST(INFO AS DOUBLE), 1.0) >= 0.8""")

# HS genome-wide significant loci (distance clumping, 500 kb)
hs = con.execute(f"SELECT SNP, CAST(CHR AS VARCHAR) CHR, POS, P FROM read_csv('{OUT}/HS.tsv.gz', delim='\\t', header=true) WHERE P < 5e-8 ORDER BY P").df()
leads = []
for r in hs.itertuples():
    if not any(l[1] == r.CHR and abs(l[2] - r.POS) < 500_000 for l in leads):
        leads.append((r.SNP, r.CHR, int(r.POS), r.P))
import pandas as pd
pd.DataFrame(leads, columns=["SNP", "CHR", "POS_b38", "P"]).to_csv(BASE / "results" / "hs_gws_loci.tsv", sep="\t", index=False)
print(f"HS GWS loci: {len(leads)}", flush=True)
# LAVA-format files: exactly one column per statistic type (SNP A1 A2 N Z), as the LAVA README requires
LAVA_IN = BASE / "data" / "lava" / "sumstats"; LAVA_IN.mkdir(parents=True, exist_ok=True)
for t in ("PSD_META_EUR", "HS", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI"):
    o = LAVA_IN / f"{t}.tsv.gz"
    if not o.exists():
        con.execute(f"COPY (SELECT SNP, A1, A2, N, Z FROM read_csv('{OUT}/{t}.tsv.gz', delim='\\t', header=true)) TO '{o}' (DELIMITER '\\t', HEADER true, COMPRESSION gzip)")
    print(f"LAVA input {t} -> {o}", flush=True)
print("SHARED_INPUTS_DONE", flush=True)
