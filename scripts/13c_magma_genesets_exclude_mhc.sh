#!/bin/bash
# Step 13c: MHC-excluded scDRS gene sets. Drops every gene overlapping the MHC (GRCh37 chr6:28.4-33.5 Mb, NCBI37.3 gene loc)
# from the MAGMA z-score table of step 13a, then re-selects the top 1,000 genes per trait with the same munge-gs settings,
# so HLA class II / antigen-presentation genes cannot drive immune cell-type associations.
# Output: data/psd_meta_magma_scdrs_noMHC.gs
set -euo pipefail
BASE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python3}
M=$BASE/data/magma
$PY - <<PYEOF
import pandas as pd
M = "$M"
loc = pd.read_csv(f"{M}/ref/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE","CHR","START","STOP","STRAND","SYMBOL"])
mhc = set(loc[(loc.CHR.astype(str) == "6") & (loc.STOP > 28_400_000) & (loc.START < 33_500_000)].SYMBOL)
z = pd.read_csv(f"{M}/out/magma_zscores.tsv", sep="\t", index_col=0)
z[~z.index.isin(mhc)].to_csv(f"{M}/out/magma_zscores_noMHC.tsv", sep="\t")
print(f"removed {z.index.isin(mhc).sum()} MHC genes of {len(z)}")
PYEOF
scdrs munge-gs --out-file $BASE/data/psd_meta_magma_scdrs_noMHC.gs --zscore-file $M/out/magma_zscores_noMHC.tsv --weight zscore --n-max 1000
awk -F'\t' 'NR>1{n=split($2,a,","); print $1, n" genes"}' $BASE/data/psd_meta_magma_scdrs_noMHC.gs
