#!/bin/bash
# Step 11a: Stream-harmonize PSD GWAS summary statistics (all GRCh38) into a common compact format.
# Output columns: chr pos a1 a2 beta se p eaf n rsid
#   a1 = effect allele after canonical ordering (a1 < a2 lexicographically; beta sign and EAF flipped if needed)
#   beta = log(OR) scale for all cohorts. Filters: autosomes 1-22, MAF >= 0.005, finite SE > 0, INFO >= 0.3 if present.
# Usage: 11a_harmonize_psd_cohorts.sh <cohort_name> [<cohort_name> ...]   (run cohorts in parallel with xargs -P)
set -euo pipefail
BASE=$(cd "$(dirname "$0")/.." && pwd)
RAW=$BASE/data/raw/sumstats
OUT=$BASE/data/processed/harmonized
mkdir -p "$OUT"

# awk function library shared by all cohorts: canonicalize alleles & print
AWK_LIB='
function emit(chr,pos,ea,oa,beta,se,p,eaf,n,rsid) {
  if (chr !~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/) return;
  sub(/^chr/,"",chr);
  if (se=="" || se=="#NA" || se=="NA" || se+0<=0) return;
  if (beta=="" || beta=="#NA" || beta=="NA") return;
  if (eaf=="" || eaf=="#NA" || eaf=="NA") return;
  eaf+=0; if (eaf<0.005 || eaf>0.995) return;
  ea=toupper(ea); oa=toupper(oa);
  if (ea !~ /^[ACGT]+$/ || oa !~ /^[ACGT]+$/ || ea==oa) return;
  if (ea < oa) { a1=ea; a2=oa; b=beta+0; f=eaf } else { a1=oa; a2=ea; b=-(beta+0); f=1-eaf }
  printf "%s\t%s\t%s\t%s\t%.6g\t%.6g\t%s\t%.5g\t%s\t%s\n", chr,pos,a1,a2,b,se,p,f,n,rsid;
}'

harmonize() {
  c=$1
  case $c in
    MVP_MultiAncestry) f=$RAW/PSD_MVP_MultiAncestry_GCST90480439.tsv.gz; mode=mvp ;;
    MVP_EUR)           f=$RAW/PSD_MVP_EUR_GCST90478766.tsv.gz;           mode=mvp ;;
    UKB_EUR)           f=$RAW/PSD_UKB_EUR_GCST90473908.tsv.gz;           mode=ukb ;;
    FinnGen_R13)       f=$RAW/PSD_FinnGen_R13_L12_PILONIDALCYST.gz;      mode=fg ;;
    GH_SAS)            f=$RAW/PSD_GH_SAS_GCST90727025.tsv;               mode=gh ;;
    *) echo "unknown cohort $c" >&2; return 1 ;;
  esac
  o=$OUT/PSD_${c}.harmonized.tsv.gz
  echo "[$(date +%T)] harmonizing $c from $f"
  if [[ $f == *.gz ]]; then reader="gzip -dc"; else reader="cat"; fi
  case $mode in
    mvp) # chr bp ea oa OR SE eaf p rsid ci_upper ci_lower alt n case_af ncase control_af nctrl r2 q i2 direction
      $reader "$f" | awk -F'\t' -v OFS='\t' "$AWK_LIB"'
        NR==1{next}
        { if ($5=="#NA"||$10=="#NA"||$11=="#NA"||$5+0<=0||$10+0<=0||$11+0<=0) next;
          beta=log($5); se=(log($10)-log($11))/(2*1.959964);
          if (se<=0) next;
          # append per-ancestry direction string to rsid field for MVP (rsid|direction|n)
          emit($1,$2,$3,$4,beta,se,$8,$7,$13,$9 "|" $21) }' ;;
    ukb) # chr bp ea oa beta se eaf p ID INFO n
      $reader "$f" | awk -F'\t' -v OFS='\t' "$AWK_LIB"'
        NR==1{next}
        { if ($10!="#NA" && $10!="" && $10+0<0.3) next; emit($1,$2,$3,$4,$5,$6,$8,$7,$11,"NA") }' ;;
    fg)  # chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
      $reader "$f" | awk -F'\t' -v OFS='\t' "$AWK_LIB"'
        NR==1{next}
        { rs=$5; if (rs=="") rs="NA"; emit($1,$2,$4,$3,$9,$10,$7,$11,"472290",rs) }' ;;
    gh)  # chr bp ea oa beta se eaf p variant_id n
      $reader "$f" | awk -F'\t' -v OFS='\t' "$AWK_LIB"'
        NR==1{next}
        { emit($1,$2,$3,$4,$5,$6,$8,$7,$10,"NA") }' ;;
  esac | { echo -e "chr\tpos\ta1\ta2\tbeta\tse\tp\teaf\tn\trsid"; cat; } | gzip -1 > "$o.tmp" && mv "$o.tmp" "$o"
  echo "[$(date +%T)] done $c -> $o ($(gzip -dc "$o" | wc -l) lines)"
}

for c in "$@"; do harmonize "$c" & done
wait
echo HARMONIZE_FINISHED "$@"
