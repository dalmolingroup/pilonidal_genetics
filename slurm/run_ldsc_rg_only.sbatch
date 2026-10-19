#!/bin/bash
#SBATCH --job-name=psd_ldsc
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/psd_ldsc_%j.out
#SBATCH --error=slurm/psd_ldsc_%j.err
# LDSC inputs (skips existing) -> munge -> h2 -> rg (PSD vs follicular-occlusion / hair / adiposity traits)
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
export DUCKDB_MEM=48GB DUCKDB_THREADS=16 LDSC_ENV=$ENVS/ldsc
$ENVS/psd/bin/python scripts/12a_prepare_ldsc_inputs.py && bash scripts/12b_run_ldsc_rg.sh || echo "LDSC FAILED"
echo "=== ldsc finished $(date) ==="
