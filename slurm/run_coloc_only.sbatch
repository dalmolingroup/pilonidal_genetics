#!/bin/bash
#SBATCH --job-name=psd_coloc
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/psd_coloc_%j.out
#SBATCH --error=slurm/psd_coloc_%j.err
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
export PATH=/mnt/beegfs/scratch/$USER/envs/rgen/bin:$PATH
Rscript scripts/14b_coloc_psd_vs_traits.R || echo "COLOC FAILED"
echo "=== coloc finished $(date) ==="
