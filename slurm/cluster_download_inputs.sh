#!/bin/bash
# Download all public inputs for Directions 2 & 3 into the repo's data/ tree (idempotent; 6 parallel streams).
set -uo pipefail
BASE=${1:-$(cd "$(dirname "$0")/.." && pwd)}
mkdir -p "$BASE/data/raw/sumstats" "$BASE/data/raw/1kg_hg38" "$BASE/data/raw/neale" "$BASE/data/ldsc_ref/eur_w_ld_chr" "$BASE/bin"
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
FG=https://storage.googleapis.com/finngen-public-data-r13/summary_stats
cat > "$BASE/data/raw/download_manifest.txt" <<U
sumstats/PSD_MVP_MultiAncestry_GCST90480439.tsv.gz $B/GCST90480001-GCST90481000/GCST90480439/GCST90480439.tsv.gz
sumstats/PSD_MVP_EUR_GCST90478766.tsv.gz $B/GCST90478001-GCST90479000/GCST90478766/GCST90478766.tsv.gz
sumstats/PSD_UKB_EUR_GCST90473908.tsv.gz $B/GCST90473001-GCST90474000/GCST90473908/GCST90473908.tsv.gz
sumstats/PSD_GH_SAS_GCST90727025.tsv $B/GCST90727001-GCST90728000/GCST90727025/GCST90727025.tsv
sumstats/PSD_FinnGen_R13_L12_PILONIDALCYST.gz $FG/finngen_R13_L12_PILONIDALCYST.gz
sumstats/HS_FinnGen_R13_L12_HIDRADENITISSUP.gz $FG/finngen_R13_L12_HIDRADENITISSUP.gz
sumstats/ACNE_FinnGen_R13_L12_ACNE.gz $FG/finngen_R13_L12_ACNE.gz
sumstats/AA_FinnGen_R13_L12_ALOPECAREATA.gz $FG/finngen_R13_L12_ALOPECAREATA.gz
sumstats/WHRadjBMI_GIANT_UKB_Pulit2019.txt.gz https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
sumstats/BMI_GIANT_UKB_Pulit2019.txt.gz https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
sumstats/HS_MVP_UKB_FinnGen_GCST90726609.h.tsv.gz $B/GCST90726001-GCST90727000/GCST90726609/harmonised/GCST90726609.h.tsv.gz
sumstats/AGA_UKB_GCST90043616.h.tsv.gz $B/GCST90043001-GCST90044000/GCST90043616/harmonised/34737426-GCST90043616-EFO_0004191.h.tsv.gz
sumstats/ACNE_TederLaving2024_GCST90245818.h.tsv.gz $B/GCST90245001-GCST90246000/GCST90245818/harmonised/GCST90245818.h.tsv.gz
1kg_hg38/all_hg38.pgen.zst https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
1kg_hg38/all_hg38.pvar.zst https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1
1kg_hg38/hg38_corrected.psam https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1
neale/2395_4.gwas.imputed_v3.male.tsv.bgz https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/2395_4.gwas.imputed_v3.male.tsv.bgz
neale/variants.tsv.bgz https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz
U
M=https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr
{ echo "ldsc_ref/eur_w_ld_chr/w_hm3.snplist $M/w_hm3.snplist"; for c in $(seq 1 22); do echo "ldsc_ref/eur_w_ld_chr/$c.l2.ldscore.gz $M/$c.l2.ldscore.gz"; echo "ldsc_ref/eur_w_ld_chr/$c.l2.M_5_50 $M/$c.l2.M_5_50"; done; } >> "$BASE/data/raw/download_manifest.txt"
dl() { rel=$1; url=$2; dest="$BASE/data/raw/$rel"; [[ $rel == ldsc_ref/* ]] && dest="$BASE/data/$rel"
  if [ -f "$dest.done" ] && [ -s "$dest" ]; then return; fi
  if curl -sSL --retry 6 --retry-delay 15 -C - -o "$dest" "$url"; then touch "$dest.done"; echo "DONE $rel $(du -h "$dest" | cut -f1)"; else echo "FAIL $rel"; fi; }
export -f dl; export BASE
xargs -P 6 -L 1 bash -c 'dl "$0" "$1"' < "$BASE/data/raw/download_manifest.txt"
# plink2 static binary
if [ ! -x "$BASE/bin/plink2" ]; then (cd "$BASE/bin" && curl -sSL -o plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20261015.zip && unzip -oq plink2.zip && ./plink2 --version); fi
echo DOWNLOADS_FINISHED
