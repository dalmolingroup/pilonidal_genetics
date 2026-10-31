#!/bin/bash
#SBATCH --job-name=psd_magma
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --output=slurm/psd_magma_%j.out
#SBATCH --error=slurm/psd_magma_%j.err
# Step 03a: MAGMA gene scores for the PSD meta-analyses and comparator traits.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
module load softwares/python/3.13.7-gnu8   # venv python links against this module's libpython3.13.so
if [ ! -x $ENVS/magma/bin/python ]; then
  python3 -m venv $ENVS/magma && $ENVS/magma/bin/pip install -q --upgrade pip && \
  $ENVS/magma/bin/pip install -q numpy==2.5.4 pandas==3.0.6 statsmodels==0.15.0 duckdb || echo "[env] magma venv install problem"
fi
$ENVS/magma/bin/python -c "import pandas, statsmodels, duckdb; print('magma env ok')"
# wait (max 60 min) for the concurrently running LDSC job to finish writing comparator inputs
for i in $(seq 1 60); do [ -s data/processed/ldsc_in/BMI_Pulit2019.tsv ] && [ -s data/processed/ldsc_in/WHRadjBMI_Pulit2019.tsv ] && break; sleep 60; done
sleep 30; ls -la data/processed/ldsc_in/*Pulit2019.tsv 2>/dev/null
PY=$ENVS/magma/bin/python bash scripts/03a_magma_gene_analysis.sh || { echo "MAGMA FAILED"; exit 1; }
echo "=== magma finished $(date) ==="
