#!/bin/bash
#SBATCH --job-name=psd_d23
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=20:00:00
#SBATCH --output=slurm/psd_d23_%j.out
#SBATCH --error=slurm/psd_d23_%j.err
# Directions 2 & 3 end-to-end on the cluster: env setup -> downloads -> harmonize -> IVW meta -> fine-mapping -> LDSC rg.
# Every stage is idempotent (skips when its outputs exist) so the job can simply be resubmitted after a failure.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD
export DUCKDB_MEM=150GB DUCKDB_THREADS=${SLURM_CPUS_PER_TASK:-32} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export PLINK2=$BASE/bin/plink2
ENVS=/mnt/beegfs/scratch/$USER/envs; mkdir -p $ENVS
echo "=== PSD Directions 2 & 3 on $(hostname) $(date) ==="

# ---- 0. environments ----
if [ ! -x $ENVS/psd/bin/python ]; then
  echo "[env] creating python venv"; ~/miniconda3/bin/python3 -m venv $ENVS/psd && $ENVS/psd/bin/pip install -q --upgrade pip && \
  $ENVS/psd/bin/pip install -q "numpy<2.3" pandas scipy matplotlib seaborn duckdb pyarrow || echo "[env] venv install problem"
fi
PY=$ENVS/psd/bin/python
if [ ! -x $ENVS/ldsc/bin/munge_sumstats.py ]; then
  echo "[env] creating ldsc conda env (python 2.7, bioconda)"; ~/miniconda3/bin/conda create -y -q -p $ENVS/ldsc -c bioconda -c conda-forge ldsc > slurm/ldsc_env_create.log 2>&1 || echo "[env] ldsc env creation FAILED (see slurm/ldsc_env_create.log)"
fi
export LDSC_ENV=$ENVS/ldsc
module load softwares/R/4.5.0-gnu11.2
export R_LIBS_USER=$ENVS/Rlib; mkdir -p $R_LIBS_USER
Rscript -e 'if(!requireNamespace("susieR",quietly=TRUE)) install.packages("susieR", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"), Ncpus=8)' > slurm/r_susie_install.log 2>&1
Rscript -e 'cat("susieR", requireNamespace("susieR",quietly=TRUE), "\n")'

# ---- 1. downloads ----
bash slurm/cluster_download_inputs.sh $BASE 2>&1 | tail -30

# ---- 2. harmonize (5 cohorts in parallel) ----
H=data/processed/harmonized; mkdir -p $H
need=""; for c in MVP_MultiAncestry MVP_EUR UKB_EUR FinnGen_R13 GH_SAS; do [ -s $H/PSD_$c.harmonized.tsv.gz ] || need="$need $c"; done
[ -n "$need" ] && bash scripts/01a_harmonize_psd_cohorts.sh $need

# ---- 3. meta-analysis ----
[ -s results/psd_meta_ALL_loci.tsv ] || $PY scripts/01b_meta_analyze_psd_gwas.py

# ---- 4. fine-mapping ----
[ -s results/psd_finemap_locus_summary.tsv ] || $PY scripts/01c_finemap_psd_loci.py

# ---- 5. LDSC inputs + rg ----
[ -s data/processed/ldsc_in/BMI_Pulit2019.tsv ] || $PY scripts/02a_prepare_ldsc_inputs.py
bash scripts/02b_run_ldsc_rg.sh
echo "=== finished $(date) ==="
