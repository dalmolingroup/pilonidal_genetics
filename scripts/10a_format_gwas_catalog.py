#!/usr/bin/env python3
"""
Step 10a: write the two PSD meta-analyses in GWAS-SSF format for GWAS Catalog deposition (run on the cluster).

Input:   data/processed/meta/psd_meta_{ALL,EUR}.tsv.gz (step 01b)
         data/processed/indel_rsids_grch38.json (step 09b; rsIDs for lead indels missing from the cohort files)
Output:  data/processed/gwas_catalog/psd_meta_{ALL,EUR}_GRCh38.tsv.gz, md5sums.txt, variant_counts.tsv

Only variants measured in >= 2 cohorts are written, the set used for discovery in the paper. Mandatory GWAS-SSF columns
come first (chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency,
p_value), then rs_id and non-standard columns described in the submission readme (n_cohorts, direction, het_q, het_p,
het_i2, n_effective). Missing values are written as #NA. Coordinates are GRCh38, 1-based.
"""
import hashlib, json, os
from pathlib import Path
import duckdb, pandas as pd

BASE = Path(__file__).resolve().parents[1]
META_DIR = BASE / "data/processed/meta"
OUT = BASE / "data/processed/gwas_catalog"
OUT.mkdir(parents=True, exist_ok=True)

# rsIDs for indels filled from Ensembl in step 09b, keyed "chr pos . a b . . ."
indels = json.load(open(BASE / "data/processed/indel_rsids_grch38.json"))
fill = pd.DataFrame([dict(chr=int(k.split()[0]), pos=int(k.split()[1]), x=k.split()[3], y=k.split()[4], rsid_fill=v)
                     for k, v in indels.items()])

con = duckdb.connect()
con.execute(f"SET threads={os.environ.get('NCPU', '4')}")
con.register("fill", fill)
counts, md5s = [], []
for meta in ("ALL", "EUR"):
    src, dst = META_DIR / f"psd_meta_{meta}.tsv.gz", OUT / f"psd_meta_{meta}_GRCh38.tsv.gz"
    con.execute(f"""
        CREATE OR REPLACE TABLE m AS
        SELECT CAST(d.chr AS INTEGER) AS chromosome, CAST(d.pos AS BIGINT) AS base_pair_location,
               a1 AS effect_allele, a2 AS other_allele,
               beta_meta AS beta, se_meta AS standard_error, eaf_meta AS effect_allele_frequency, p_meta AS p_value,
               COALESCE(NULLIF(rsid, ''), f.rsid_fill) AS rs_id,
               k AS n_cohorts, direction, q_het AS het_q, p_het AS het_p, i2 AS het_i2, n_eff AS n_effective
        FROM read_csv('{src}', delim='\t', header=true, all_varchar=false,
                      types={{'chr': 'VARCHAR', 'rsid': 'VARCHAR', 'direction': 'VARCHAR', 'mvp_direction': 'VARCHAR'}}) d
        LEFT JOIN fill f ON CAST(d.chr AS INTEGER) = f.chr AND d.pos = f.pos
             AND ((d.a1 = f.x AND d.a2 = f.y) OR (d.a1 = f.y AND d.a2 = f.x))
        WHERE k >= 2
    """)
    c = con.execute("""SELECT count(*) n, count(*) FILTER (WHERE effect_allele_frequency IS NULL) eaf_na,
                              count(*) FILTER (WHERE rs_id IS NULL) rsid_na, min(p_value) min_p,
                              count(*) FILTER (WHERE p_value = 0) p_zero, count(*) FILTER (WHERE p_value < 5e-8) gws
                       FROM m""").df().assign(meta=f"META_{meta}")
    counts.append(c)
    print(c.to_string(index=False), flush=True)
    con.execute(f"""COPY (SELECT * FROM m ORDER BY chromosome, base_pair_location)
                    TO '{dst}' (FORMAT CSV, DELIMITER '\t', HEADER, NULLSTR '#NA', COMPRESSION GZIP)""")
    h = hashlib.md5()
    with open(dst, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    md5s.append(f"{h.hexdigest()}  {dst.name}")
    print(f"wrote {dst} ({dst.stat().st_size / 1e6:.0f} MB)", flush=True)

(OUT / "md5sums.txt").write_text("\n".join(md5s) + "\n")
pd.concat(counts).to_csv(OUT / "variant_counts.tsv", sep="\t", index=False)
print("\n".join(md5s))
