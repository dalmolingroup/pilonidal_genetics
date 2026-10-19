#!/usr/bin/env python3
"""
Step 12a: Build munge-ready summary statistics (HapMap3 SNPs only) for LDSC genetic correlation analyses.
Produces data/processed/ldsc_in/<trait>.tsv with columns SNP A1 A2 BETA SE P N  (A1 = effect allele).

PSD sets: META_EUR, META_ALL, MVP_EUR, MVP_MultiAncestry, UKB_EUR, FinnGen_R13
Traits:   HS (GCST90726609 MVP+UKB+FinnGen; FinnGen R13), acne (GCST90245818 Teder-Laving 2024; FinnGen R13),
          male-pattern baldness (GCST90043616 UKB; Neale UKB 2395_4 pattern 4), alopecia areata (FinnGen R13),
          WHRadjBMI and BMI (Pulit 2019 GIANT+UKB).
"""
import os
from pathlib import Path
import duckdb

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "data" / "raw"
S = RAW / "sumstats"
HARM = BASE / "data" / "processed" / "harmonized"
META = BASE / "data" / "processed" / "meta"
OUT = BASE / "data" / "processed" / "ldsc_in"
OUT.mkdir(parents=True, exist_ok=True)
HM3 = BASE / "data" / "ldsc_ref" / "eur_w_ld_chr" / "w_hm3.snplist"

con = duckdb.connect()
con.execute(f"SET memory_limit='{os.environ.get('DUCKDB_MEM','12GB')}'; SET threads={os.environ.get('DUCKDB_THREADS','12')};")
con.execute(f"SET temp_directory='{OUT}/duckdb_tmp';")
con.execute(f"CREATE TABLE hm3 AS SELECT SNP FROM read_csv('{HM3}', delim='\\t', header=true)")
print("HM3 SNPs:", con.execute("SELECT COUNT(*) FROM hm3").fetchone()[0], flush=True)

# rsID map for GRCh38 keys (chr,pos,a1,a2 canonical) from MVP + FinnGen harmonized files, restricted to HM3
con.execute(f"""
CREATE TABLE rsmap AS
SELECT chr, CAST(pos AS BIGINT) AS pos, a1, a2, ANY_VALUE(rs) AS SNP FROM (
  SELECT chr, pos, a1, a2, split_part(rsid,'|',1) AS rs FROM read_csv('{HARM}/PSD_MVP_MultiAncestry.harmonized.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
  UNION ALL
  SELECT chr, pos, a1, a2, unnest(string_split(rsid, ',')) AS rs FROM read_csv('{HARM}/PSD_FinnGen_R13.harmonized.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
) WHERE rs IN (SELECT SNP FROM hm3) GROUP BY chr, pos, a1, a2
""")
print("rsmap HM3 keys:", con.execute("SELECT COUNT(*) FROM rsmap").fetchone()[0], flush=True)


def write(name, sql):
    out = OUT / f"{name}.tsv"
    if out.exists() and out.stat().st_size > 1_000_000:
        print(f"  {name}: exists, skipping", flush=True); return
    con.execute(f"COPY ({sql}) TO '{out}' (DELIMITER '\\t', HEADER true)")
    n = con.execute(f"SELECT COUNT(*) FROM read_csv('{out}', delim='\\t', header=true)").fetchone()[0]
    print(f"  {name}: {n:,} HM3 SNPs -> {out}", flush=True)


# ---- PSD meta-analyses (rsid from meta table; fill missing via rsmap) ----
for m, n_total in (("EUR", 1433 + 448820 + 1466 + 456974 + 3872 + 468418), ("ALL", 2073 + 629531 + 1466 + 456974 + 3872 + 468418 + 457 + 43569)):
    write(f"PSD_META_{m}", f"""
      SELECT COALESCE(r.SNP, t.rsid) AS SNP, t.a1 AS A1, t.a2 AS A2, CAST(t.beta_meta AS DOUBLE) AS BETA, CAST(t.se_meta AS DOUBLE) AS SE,
             CAST(t.p_meta AS DOUBLE) AS P, CAST(ROUND({n_total} * CAST(t.n_eff AS DOUBLE) / MAX(CAST(t.n_eff AS DOUBLE)) OVER ()) AS BIGINT) AS N
      FROM read_csv('{META}/psd_meta_{m}.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) t
      LEFT JOIN rsmap r ON r.chr=t.chr AND r.pos=CAST(t.pos AS BIGINT) AND r.a1=t.a1 AND r.a2=t.a2
      WHERE CAST(t.k AS INT) >= 2 AND COALESCE(r.SNP, t.rsid) IN (SELECT SNP FROM hm3)
    """)

# ---- individual PSD cohorts ----
for c, n_total in (("MVP_EUR", 450253), ("MVP_MultiAncestry", 631604), ("UKB_EUR", 458440), ("FinnGen_R13", 472290)):
    write(f"PSD_{c}", f"""
      SELECT r.SNP, t.a1 AS A1, t.a2 AS A2, CAST(t.beta AS DOUBLE) AS BETA, CAST(t.se AS DOUBLE) AS SE, TRY_CAST(t.p AS DOUBLE) AS P,
             CAST(COALESCE(TRY_CAST(t.n AS DOUBLE), {n_total}) AS BIGINT) AS N
      FROM read_csv('{HARM}/PSD_{c}.harmonized.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) t
      JOIN rsmap r ON r.chr=t.chr AND r.pos=CAST(t.pos AS BIGINT) AND r.a1=t.a1 AND r.a2=t.a2
      WHERE CAST(t.eaf AS DOUBLE) BETWEEN 0.01 AND 0.99
    """)

# ---- Hidradenitis suppurativa: GCST90726609 (MVP+UKB+FinnGen, 3,941 cases / 1,435,603 controls) ----
write("HS_GCST90726609", f"""
  SELECT rsid AS SNP, effect_allele AS A1, other_allele AS A2, CAST(beta AS DOUBLE) AS BETA, CAST(standard_error AS DOUBLE) AS SE,
         CAST(p_value AS DOUBLE) AS P, 1439544 AS N
  FROM read_csv('{S}/HS_MVP_UKB_FinnGen_GCST90726609.h.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
  WHERE rsid IN (SELECT SNP FROM hm3) AND TRY_CAST(standard_error AS DOUBLE) > 0 AND CAST(effect_allele_frequency AS DOUBLE) BETWEEN 0.01 AND 0.99
""")

# ---- Acne: GCST90245818 (Teder-Laving 2024; 34,422 cases / 364,991 controls) ----
write("ACNE_GCST90245818", f"""
  SELECT rsid AS SNP, effect_allele AS A1, other_allele AS A2, CAST(beta AS DOUBLE) AS BETA, CAST(standard_error AS DOUBLE) AS SE,
         CAST(p_value AS DOUBLE) AS P, 399413 AS N
  FROM read_csv('{S}/ACNE_TederLaving2024_GCST90245818.h.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
  WHERE rsid IN (SELECT SNP FROM hm3) AND TRY_CAST(standard_error AS DOUBLE) > 0 AND TRY_CAST(effect_allele_frequency AS DOUBLE) BETWEEN 0.01 AND 0.99
""")

# ---- Male-pattern baldness: GCST90043616 (UKB; 66,172 cases / 140,864 controls) ----
write("MPB_GCST90043616", f"""
  SELECT hm_rsid AS SNP, hm_effect_allele AS A1, hm_other_allele AS A2, CAST(hm_beta AS DOUBLE) AS BETA, CAST(standard_error AS DOUBLE) AS SE,
         CAST(p_value AS DOUBLE) AS P, CAST(COALESCE(TRY_CAST(n AS DOUBLE), 207036) AS BIGINT) AS N
  FROM read_csv('{S}/AGA_UKB_GCST90043616.h.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
  WHERE hm_rsid IN (SELECT SNP FROM hm3) AND TRY_CAST(standard_error AS DOUBLE) > 0 AND TRY_CAST(hm_effect_allele_frequency AS DOUBLE) BETWEEN 0.01 AND 0.99
""")

# ---- Neale UKB 2395_4 Hair/balding pattern 4 (males; 30,225 cases / 135,424 controls) ----
write("BALDING4_Neale_UKB_male", f"""
  SELECT v.rsid AS SNP, v.alt AS A1, v.ref AS A2, CAST(g.beta AS DOUBLE) AS BETA, CAST(g.se AS DOUBLE) AS SE, CAST(g.pval AS DOUBLE) AS P,
         CAST(g.n_complete_samples AS BIGINT) AS N
  FROM read_csv('{RAW}/neale/2395_4.gwas.imputed_v3.male.tsv.bgz', compression='gzip', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) g
  JOIN read_csv('{RAW}/neale/variants.tsv.bgz', compression='gzip', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) v ON v.variant = g.variant
  WHERE v.rsid IN (SELECT SNP FROM hm3) AND g.low_confidence_variant = 'false' AND CAST(g.minor_AF AS DOUBLE) >= 0.01 AND TRY_CAST(g.se AS DOUBLE) > 0
""")

# ---- FinnGen R13 skin endpoints ----
for name, f, n in (("HS_FinnGen_R13", "HS_FinnGen_R13_L12_HIDRADENITISSUP.gz", 1453 + 474896),
                   ("ACNE_FinnGen_R13", "ACNE_FinnGen_R13_L12_ACNE.gz", 4812 + 474896),
                   ("AA_FinnGen_R13", "AA_FinnGen_R13_L12_ALOPECAREATA.gz", 994 + 474896)):
    write(name, f"""
      SELECT rs AS SNP, alt AS A1, ref AS A2, CAST(beta AS DOUBLE) AS BETA, CAST(sebeta AS DOUBLE) AS SE, CAST(pval AS DOUBLE) AS P, {n} AS N
      FROM (SELECT *, unnest(string_split(rsids, ',')) AS rs FROM read_csv('{S}/{f}', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000))
      WHERE rs IN (SELECT SNP FROM hm3) AND TRY_CAST(sebeta AS DOUBLE) > 0 AND CAST(af_alt AS DOUBLE) BETWEEN 0.01 AND 0.99
    """)

# ---- GIANT+UKB Pulit 2019 ----
for name, f in (("WHRadjBMI_Pulit2019", "WHRadjBMI_GIANT_UKB_Pulit2019.txt.gz"), ("BMI_Pulit2019", "BMI_GIANT_UKB_Pulit2019.txt.gz")):
    write(name, f"""
      SELECT split_part(SNP, ':', 1) AS SNP, Tested_Allele AS A1, Other_Allele AS A2, TRY_CAST(BETA AS DOUBLE) AS BETA, TRY_CAST(SE AS DOUBLE) AS SE,
             TRY_CAST(P AS DOUBLE) AS P, CAST(TRY_CAST(N AS DOUBLE) AS BIGINT) AS N
      FROM read_csv('{S}/{f}', delim=' ', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000)
      WHERE split_part(SNP, ':', 1) IN (SELECT SNP FROM hm3) AND TRY_CAST(Freq_Tested_Allele AS DOUBLE) BETWEEN 0.01 AND 0.99 AND COALESCE(TRY_CAST(INFO AS DOUBLE), 1.0) >= 0.8 AND TRY_CAST(SE AS DOUBLE) > 0 AND TRY_CAST(P AS DOUBLE) IS NOT NULL AND TRY_CAST(N AS DOUBLE) > 0 AND TRY_CAST(BETA AS DOUBLE) IS NOT NULL
    """)
print("LDSC_INPUTS_DONE", flush=True)
