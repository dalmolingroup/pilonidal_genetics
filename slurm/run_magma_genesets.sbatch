#!/bin/bash
#SBATCH --job-name=psd_magma_gsa
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/psd_magma_gsa_%j.out
#SBATCH --error=slurm/psd_magma_gsa_%j.err
# Step 06a: MHC-free MAGMA gene-set analyses (GO BP + Reactome) for PSD and comparators.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
module load softwares/python/3.13.7-gnu8   # magma venv python links against this module's libpython
PY=$ENVS/magma/bin/python bash scripts/06a_magma_genesets.sh || echo "MAGMA GENESETS FAILED"
echo "=== magma gene-set analysis finished $(date) ==="
