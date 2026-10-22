#!/bin/bash
#SBATCH --job-name=psd_rerun
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/psd_rerun_%j.out
#SBATCH --error=slurm/psd_rerun_%j.err
# Re-run after the per-cohort duplicate-key fix in 01b: meta -> fine-mapping (with SuSiE from conda-forge) -> LDSC inputs/rg.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD; ENVS=/mnt/beegfs/scratch/$USER/envs
export DUCKDB_MEM=150GB DUCKDB_THREADS=${SLURM_CPUS_PER_TASK:-32} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export PLINK2=$BASE/bin/plink2 LDSC_ENV=$ENVS/ldsc
PY=$ENVS/psd/bin/python
if [ ! -x $ENVS/rsusie/bin/Rscript ]; then
  echo "[env] creating conda env with r-susier"; ~/miniconda3/bin/conda create -y -q -p $ENVS/rsusie -c conda-forge r-base=4.3 r-susier r-data.table > slurm/rsusie_env_create.log 2>&1 || echo "[env] r-susier env creation FAILED"
fi
export PATH=$ENVS/rsusie/bin:$PATH
Rscript -e 'cat("susieR", requireNamespace("susieR",quietly=TRUE), "\n")'
# Replace the 2021-preliminary (b38a) pvar with the 2022-08-04 (b38b) pvar that matches the downloaded pgen
KG=data/raw/1kg_hg38
if [ ! -f $KG/.pvar_b38b ]; then rm -f $KG/all_hg38.pvar.zst $KG/all_hg38.pvar.zst.done; bash slurm/cluster_download_inputs.sh $BASE 2>&1 | grep -E "pvar|FAIL"; touch $KG/.pvar_b38b; fi
echo "pvar variants: $(zstd -dc $KG/all_hg38.pvar.zst | grep -vc '^#')  (pgen expects 75193455)"
mkdir -p results/run1_before_dedup && cp -f results/psd_meta_* results/psd_finemap_* results/psd_published_loci_replication.tsv results/run1_before_dedup/ 2>/dev/null || true
rm -f results/psd_meta_ALL_loci.tsv results/psd_meta_EUR_loci.tsv results/psd_finemap_locus_summary.tsv
rm -rf data/processed/ldsc_in data/processed/ldsc_munged results/ldsc_logs
$PY scripts/01b_meta_analyze_psd_gwas.py || { echo "META FAILED"; exit 1; }
$PY scripts/01c_finemap_psd_loci.py || echo "FINEMAP FAILED (continuing to LDSC)"
$PY scripts/02a_prepare_ldsc_inputs.py && bash scripts/02b_run_ldsc_rg.sh || echo "LDSC FAILED"
echo "=== rerun finished $(date) ==="
