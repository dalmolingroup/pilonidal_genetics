#!/bin/bash
#SBATCH --job-name=psd_susie_ukb
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/psd_susie_ukb_%j.out
#SBATCH --error=slurm/psd_susie_ukb_%j.err
# SuSiE-RSS sensitivity analysis (01f) with UK Biobank LD; downloads ~25-50 GB of LD files into data/raw/ukbb_ld on first run.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD; ENVS=/mnt/beegfs/scratch/$USER/envs
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export PATH=$ENVS/rsusie/bin:$PATH
Rscript -e 'cat("susieR", as.character(packageVersion("susieR")), "\n")'
df -h $BASE | tail -1
$ENVS/psd/bin/python scripts/01f_finemap_susie_ukbld.py
echo "=== susie ukb-ld finished $(date) ==="
