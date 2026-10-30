#!/bin/bash
# Step 02b: munge all LDSC inputs, estimate univariate h2 for PSD sets, and genetic correlations (rg)
# between PSD and follicular-occlusion / hair / adiposity traits, plus cross-cohort PSD rg among independent cohorts.
set -euo pipefail
BASE=$(cd "$(dirname "$0")/.." && pwd)
LDSC_ENV=${LDSC_ENV:-/mnt/beegfs/scratch/$USER/envs/ldsc}
PY=$LDSC_ENV/bin/python
REF=$BASE/data/ldsc_ref/eur_w_ld_chr
IN=$BASE/data/processed/ldsc_in
MUNGED=$BASE/data/processed/ldsc_munged
LOG=$BASE/results/ldsc_logs
mkdir -p "$MUNGED" "$LOG"

munge() {
  t=$1
  if [ -s "$MUNGED/$t.sumstats.gz" ]; then return; fi
  $PY $LDSC_ENV/bin/munge_sumstats.py --sumstats "$IN/$t.tsv" --out "$MUNGED/$t" --merge-alleles "$REF/w_hm3.snplist" \
      --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-col N --chunksize 500000 > "$LOG/munge_$t.log" 2>&1 || { echo "munge FAILED $t"; tail -5 "$LOG/munge_$t.log"; }
}
export -f munge; export PY LDSC_ENV REF IN MUNGED LOG
ls "$IN"/*.tsv | xargs -n1 basename | sed 's/\.tsv$//' | xargs -P 6 -I{} bash -c 'munge {}'
echo "munged: $(ls $MUNGED/*.sumstats.gz | wc -l)"

PSD_SETS="PSD_META_EUR PSD_META_ALL PSD_MVP_EUR PSD_MVP_MultiAncestry PSD_UKB_EUR PSD_FinnGen_R13"
TRAITS="HS_GCST90726609 HS_FinnGen_R13 ACNE_GCST90245818 ACNE_FinnGen_R13 MPB_GCST90043616 BALDING4_Neale_UKB_male AA_FinnGen_R13 WHRadjBMI_Pulit2019 BMI_Pulit2019"

# Univariate h2 (observed scale) for each PSD set and trait
for t in $PSD_SETS $TRAITS; do
  [ -s "$LOG/h2_$t.log" ] && continue
  $PY $LDSC_ENV/bin/ldsc.py --h2 "$MUNGED/$t.sumstats.gz" --ref-ld-chr "$REF/" --w-ld-chr "$REF/" --out "$LOG/h2_$t" > /dev/null 2>&1 || echo "h2 failed $t"
done

# rg: each PSD set vs all traits (one ldsc call per PSD set)
for p in $PSD_SETS; do
  [ -s "$LOG/rg_${p}_vs_traits.log" ] && continue
  list="$MUNGED/$p.sumstats.gz"; for t in $TRAITS; do list="$list,$MUNGED/$t.sumstats.gz"; done
  $PY $LDSC_ENV/bin/ldsc.py --rg "$list" --ref-ld-chr "$REF/" --w-ld-chr "$REF/" --out "$LOG/rg_${p}_vs_traits" > /dev/null 2>&1 || echo "rg failed $p"
done

# Cross-cohort PSD rg among independent cohorts (MVP_EUR, UKB_EUR, FinnGen_R13)
[ -s "$LOG/rg_psd_cross_cohort.log" ] || $PY $LDSC_ENV/bin/ldsc.py --rg "$MUNGED/PSD_MVP_EUR.sumstats.gz,$MUNGED/PSD_UKB_EUR.sumstats.gz,$MUNGED/PSD_FinnGen_R13.sumstats.gz" \
    --ref-ld-chr "$REF/" --w-ld-chr "$REF/" --out "$LOG/rg_psd_cross_cohort" > /dev/null 2>&1 || echo "cross-cohort rg failed"
[ -s "$LOG/rg_psd_ukb_finngen.log" ] || $PY $LDSC_ENV/bin/ldsc.py --rg "$MUNGED/PSD_UKB_EUR.sumstats.gz,$MUNGED/PSD_FinnGen_R13.sumstats.gz" \
    --ref-ld-chr "$REF/" --w-ld-chr "$REF/" --out "$LOG/rg_psd_ukb_finngen" > /dev/null 2>&1 || true
echo LDSC_RG_DONE
