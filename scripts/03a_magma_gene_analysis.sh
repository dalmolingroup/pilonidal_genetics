#!/bin/bash
# Step 03a: MAGMA v1.10 gene-level association (SNP-wise mean model; 1000G EUR LD reference; +10 kb / -1.5 kb windows)
# for the PSD meta-analyses and comparator traits, then scDRS gene sets (top 1,000 genes by MAGMA z, z-score weights).
# PSD inputs use all rsID-mapped meta-analysis variants (k >= 2); comparators use their HapMap3 LDSC inputs.
set -euo pipefail
BASE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python3}
M=$BASE/data/magma; mkdir -p $M/ref $M/in $M/out
cd $M/ref
S=https://vu.data.surfsara.nl/index.php/s
[ -x magma ] || { curl -sSL -o magma.zip $S/lxDgt2dNdNr6DYt/download && unzip -oq magma.zip && chmod +x magma; }
[ -s g1000_eur.bed ] || { curl -sSL -o g1000_eur.zip $S/VZNByNwpD8qqINe/download && unzip -oq g1000_eur.zip; }
[ -s NCBI37.3.gene.loc ] || { curl -sSL -o geneloc37.zip $S/Pj2orwuF2JYyKxq/download && unzip -oq geneloc37.zip; }
ls -la $M/ref | head; ./magma --version || true
[ -s $M/ref/annot.genes.annot ] || ./magma --annotate window=10,1.5 --snp-loc g1000_eur.bim --gene-loc NCBI37.3.gene.loc --out $M/ref/annot > /dev/null

# ---- build SNP P N inputs ----
$PY - <<PYEOF
import duckdb, os
con = duckdb.connect(); con.execute("SET threads=8; SET memory_limit='40GB';")
B = "$BASE"; M = "$M"
opts = "delim='\\t', header=true, all_varchar=true, quote='', escape='', strict_mode=false, ignore_errors=true"
for m in ("EUR", "ALL"):
    con.execute(f"""COPY (SELECT split_part(rsid,'|',1) AS SNP, MIN(CAST(p_meta AS DOUBLE)) AS P, CAST(ROUND(MAX(CAST(n_eff AS DOUBLE))) AS BIGINT) AS N
        FROM read_csv('{B}/data/processed/meta/psd_meta_{m}.tsv.gz', {opts})
        WHERE CAST(k AS INT) >= 2 AND rsid LIKE 'rs%' GROUP BY 1) TO '{M}/in/PSD_META_{m}.txt' (DELIMITER ' ', HEADER true)""")
for t in ("HS_GCST90726609", "MPB_GCST90043616", "ACNE_GCST90245818", "WHRadjBMI_Pulit2019", "BMI_Pulit2019"):
    f = f"{B}/data/processed/ldsc_in/{t}.tsv"
    if os.path.exists(f):
        con.execute(f"""COPY (SELECT SNP, MIN(CAST(P AS DOUBLE)) P, MAX(CAST(N AS BIGINT)) N FROM read_csv('{f}', {opts}) GROUP BY 1)
                        TO '{M}/in/{t}.txt' (DELIMITER ' ', HEADER true)""")
print("MAGMA inputs written")
PYEOF

# ---- gene analysis (parallel over traits) ----
run_magma() { t=$1; [ -s $M/out/$t.genes.out ] && return; $M/ref/magma --bfile $M/ref/g1000_eur --pval $M/in/$t.txt use=SNP,P ncol=N \
    --gene-annot $M/ref/annot.genes.annot --out $M/out/$t > $M/out/$t.stdout 2>&1 || echo "MAGMA failed $t"; }
export -f run_magma; export M
ls $M/in/*.txt | xargs -n1 basename | sed 's/\.txt$//' | xargs -P 6 -I{} bash -c 'run_magma {}'

# ---- z-score table (gene symbols) + scDRS gene sets ----
$PY - <<PYEOF
import pandas as pd, glob, os
M = "$M"; B = "$BASE"
loc = pd.read_csv(f"{M}/ref/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE","CHR","START","STOP","STRAND","SYMBOL"])
z = None
for f in sorted(glob.glob(f"{M}/out/*.genes.out")):
    t = os.path.basename(f).replace(".genes.out", "")
    d = pd.read_csv(f, sep=r"\s+")[["GENE", "ZSTAT", "P"]].merge(loc[["GENE", "SYMBOL"]], on="GENE")
    d.to_csv(f"{B}/results/magma_genes_{t}.tsv", sep="\t", index=False)
    n_sig = int((d.P < 0.05 / len(d)).sum())
    print(f"{t}: {len(d):,} genes, Bonferroni-significant {n_sig}; top: {', '.join(d.sort_values('P').SYMBOL.head(12))}")
    s = d.drop_duplicates("SYMBOL").set_index("SYMBOL")["ZSTAT"].rename(t)
    z = s.to_frame() if z is None else z.join(s, how="outer")
z = z.fillna(0.0); z.index.name = "GENE"
z.to_csv(f"{M}/out/magma_zscores.tsv", sep="\t")
PYEOF
scdrs munge-gs --out-file $BASE/data/psd_meta_magma_scdrs.gs --zscore-file $M/out/magma_zscores.tsv --weight zscore --n-max 1000
awk -F'\t' 'NR>1{n=split($2,a,","); print $1, n" genes"}' $BASE/data/psd_meta_magma_scdrs.gs
echo MAGMA_DONE
