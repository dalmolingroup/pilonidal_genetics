#!/bin/bash
#SBATCH --job-name=psd_shared
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/psd_shared_%j.out
#SBATCH --error=slurm/psd_shared_%j.err
# Direction 3 extension: PSD-HS-acne-hair-adiposity shared genetics via coloc (04b), LAVA local rg (04c) and Genomic SEM (04d).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD; ENVS=/mnt/beegfs/scratch/$USER/envs
export DUCKDB_MEM=150GB DUCKDB_THREADS=32 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
LDSC_PY=$ENVS/ldsc/bin/python; LDSC=$ENVS/ldsc/bin/ldsc.py; REF=$BASE/data/ldsc_ref/eur_w_ld_chr/
echo "=== shared genetics job on $(hostname) $(date) ==="

# ---- 0. R environment (conda-forge R 4.3 + compilers; coloc/GenomicSEM/LAVA from CRAN/GitHub) ----
RE=$ENVS/rgen
if [ ! -x $RE/bin/Rscript ]; then
  ~/miniconda3/bin/conda create -y -q -p $RE -c conda-forge r-base=4.3 r-data.table r-lavaan r-matrix r-remotes r-rcpp r-dplyr r-stringr \
     r-plyr r-foreach r-doparallel r-iterators r-readr r-e1071 r-gdata r-r.utils r-susier r-mvtnorm r-ggplot2 r-viridis r-splitstackshape \
     compilers make > slurm/rgen_env_create.log 2>&1 || {
     echo "[env] full rgen creation failed; retrying with a minimal set (rest from CRAN)"; rm -rf $RE
     ~/miniconda3/bin/conda create -y -q -p $RE -c conda-forge r-base=4.3 r-data.table r-remotes r-lavaan r-rcpp compilers make >> slurm/rgen_env_create.log 2>&1 || echo "[env] rgen creation FAILED"; }
fi
Rscript_ok() { $RE/bin/Rscript -e 'r <- "https://cloud.r-project.org"; need <- c("data.table","lavaan","Matrix","remotes","Rcpp","dplyr","stringr","plyr","foreach","doParallel","iterators","readr","e1071","gdata","R.utils","susieR","mvtnorm","ggplot2","viridis","splitstackshape"); miss <- need[!sapply(need, requireNamespace, quietly=TRUE)]; if (length(miss)) install.packages(miss, repos=r, Ncpus=8)' >> slurm/rgen_pkgs.log 2>&1; }
Rscript_ok
export PATH=$RE/bin:$PATH
Rscript -e 'r <- "https://cloud.r-project.org"; for (p in c("coloc","mgsub","simsalapar","splitstackshape")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos=r, Ncpus=8);
  if (!requireNamespace("GenomicSEM", quietly=TRUE)) remotes::install_github("GenomicSEM/GenomicSEM", upgrade="never", dependencies=TRUE);
  if (!requireNamespace("LAVA", quietly=TRUE)) remotes::install_github("josefin-werme/LAVA", upgrade="never", dependencies=TRUE)' > slurm/rgen_pkgs.log 2>&1
Rscript -e 'for (p in c("coloc","GenomicSEM","LAVA","data.table")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'

# ---- 1. inputs ----
$ENVS/psd/bin/python scripts/04a_prepare_shared_genetics_inputs.py || { echo "INPUTS FAILED"; exit 1; }

# ---- 2. LAVA reference (UKB v1.1) + locus file ----
L=$BASE/data/lava; mkdir -p $L/ref
[ -s $L/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile ] || curl -sSL -o $L/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile \
   https://raw.githubusercontent.com/josefin-werme/LAVA/main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
if [ ! -f $L/ref/.complete ]; then
  ( cd $L/ref
    printf "%s\n" "7NBVIvtPRdu7Qhz chr1-2" "fy6ITboMojHrQXr chr3-4" "mRz31q0lq7KMcuI chr5-6" "89RXxPN2BlOqxLs chr7-9" "3dU1L2Hap43xuCs chr10-12" "DQqJ2Sqr49RP4xe chr13-16" "U8eg5XfTr8qPzPp chr17-23" | \
      xargs -P 7 -L 1 bash -c 'curl -sSL --retry 5 -o lava-ukb-v1.1_$1.zip https://vu.data.surfsara.nl/index.php/s/$0/download && unzip -oq lava-ukb-v1.1_$1.zip && rm -f lava-ukb-v1.1_$1.zip && echo "LAVA ref $1 ok"'
    find . -maxdepth 2 -type f | head -20; touch .complete )
fi
# flatten if the zips extracted into sub-folders
find $L/ref -mindepth 2 -type f -name 'lava-ukb-v1.1*' -exec mv -n {} $L/ref/ \; 2>/dev/null; ls $L/ref | head -12

# ---- 3. sample overlap matrix from LDSC cross-trait intercepts + LAVA info file ----
M=$BASE/data/processed/ldsc_munged; O=$BASE/results/ldsc_overlap; mkdir -p $O
declare -A F=( [PSD]=PSD_META_EUR [HS]=HS_GCST90726609 [ACNE_FG]=ACNE_FinnGen_R13 [ACNE_TL]=ACNE_GCST90245818 [BALDP1]=MPB_GCST90043616 [BMI]=BMI_Pulit2019 [WHRadjBMI]=WHRadjBMI_Pulit2019 )
ORDER="PSD HS ACNE_FG ACNE_TL BALDP1 BMI WHRadjBMI"
for t in $ORDER; do
  [ -s $O/$t.log ] && continue
  list="$M/${F[$t]}.sumstats.gz"; for u in $ORDER; do list="$list,$M/${F[$u]}.sumstats.gz"; done
  $LDSC_PY $LDSC --rg "$list" --ref-ld-chr $REF --w-ld-chr $REF --out $O/$t > /dev/null 2>&1 &
done; wait
$ENVS/psd/bin/python - <<'PYEOF'
import re, numpy as np, pandas as pd
O = "results/ldsc_overlap"; L = "data/lava"
F = dict(PSD="PSD_META_EUR", HS="HS_GCST90726609", ACNE_FG="ACNE_FinnGen_R13", ACNE_TL="ACNE_GCST90245818", BALDP1="MPB_GCST90043616", BMI="BMI_Pulit2019", WHRadjBMI="WHRadjBMI_Pulit2019")
inv = {v: k for k, v in F.items()}; names = list(F)
mat = pd.DataFrame(np.nan, index=names, columns=names)
for t in names:
    txt = open(f"{O}/{t}.log").read()
    blk = re.search(r"Summary of Genetic Correlation Results\n(.*?)\n\n", txt, re.S).group(1).strip().split("\n")
    hdr = blk[0].split()
    for line in blk[1:]:
        d = dict(zip(hdr, line.split()))
        p1 = inv[d["p1"].split("/")[-1].replace(".sumstats.gz", "")]; p2 = inv[d["p2"].split("/")[-1].replace(".sumstats.gz", "")]
        g = d.get("gcov_int", "NA")
        if g not in ("NA", "nan"): mat.loc[p1, p2] = float(g)
        if p1 == p2 and d.get("h2_int", "NA") not in ("NA", "nan"): mat.loc[p1, p1] = float(d["h2_int"])
mat = (mat + mat.T) / 2
d = np.sqrt(np.diag(mat.values)); cor = mat.values / np.outer(d, d)
cor = pd.DataFrame(np.round(cor, 5), index=names, columns=names)
cor.to_csv(f"{L}/sample.overlap.txt", sep=" ")
print(cor)
info = pd.DataFrame({"phenotype": names,
    "cases":    [6771, 3941, 4812, 34422, 66172, np.nan, np.nan],
    "controls": [1374229, 1435603, 474896, 364991, 140864, np.nan, np.nan],
    "prevalence": [0.007, 0.01, 0.01, 34422/399413, 66172/207036, np.nan, np.nan],
    "filename": [f"data/lava/sumstats/{x}.tsv.gz" for x in ["PSD_META_EUR", "HS", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI"]]})
info.to_csv(f"{L}/input.info.txt", sep="\t", index=False, na_rep="NA")
print(info)
PYEOF

# ---- 4-6. analyses ----
[ -s results/coloc_psd_vs_traits.tsv ] || Rscript scripts/04b_coloc_psd_vs_traits.R || echo "COLOC FAILED"
[ -s results/genomicsem_models.tsv ] || Rscript scripts/04d_genomic_sem.R || echo "GENOMICSEM FAILED"
[ -s results/lava_bivariate_psd.tsv ] || Rscript scripts/04c_lava_local_rg.R || echo "LAVA FAILED"
echo "=== shared genetics finished $(date) ==="
