#!/bin/bash
#SBATCH --job-name=psd_ldsc_sldsc
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/psd_ldsc_sldsc_%j.out
#SBATCH --error=slurm/psd_ldsc_sldsc_%j.err

set -euo pipefail

cd /scratch/global/$USER/pilonidal_genetics

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "=== Starting PSD Univariate LDSC, Bivariate rg, and Stratified LDSC (LDSC-SEG) on $(hostname) ==="
date

/opt/npad/shared/conda/forge2023/bin/python3 scripts/10_ldsc_heritability_and_celltype_partitioning.py

echo "=== Completed PSD LDSC & S-LDSC ==="
date
