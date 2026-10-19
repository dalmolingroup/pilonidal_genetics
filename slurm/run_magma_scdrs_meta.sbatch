#!/bin/bash
#SBATCH --job-name=psd_magma_scdrs
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --output=slurm/psd_magma_scdrs_%j.out
#SBATCH --error=slurm/psd_magma_scdrs_%j.err
# Direction 1 v2: MAGMA gene scores from the PSD meta-analysis -> scDRS across 7 skin / hair-follicle / HS atlases.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2   # must match the values script 09 sets at import
module load softwares/python/3.13.7-gnu8   # venv python links against this module's libpython3.13.so
if [ ! -x $ENVS/scdrs/bin/scdrs ]; then
  python3 -m venv $ENVS/scdrs && $ENVS/scdrs/bin/pip install -q --upgrade pip && \
  $ENVS/scdrs/bin/pip install -q numpy==2.5.4 pandas==3.0.6 scipy==1.18.1 scanpy==1.12.4 anndata==0.13.4 scdrs==1.0.2 statsmodels==0.15.0 \
      leidenalg==0.12.0 igraph==1.0.0 h5py==3.16.0 numba==0.68.0 scikit-learn==1.9.1 matplotlib seaborn duckdb || echo "[env] scdrs venv install problem"
fi
export PATH=$ENVS/scdrs/bin:$PATH
$ENVS/scdrs/bin/python -c "import scdrs, scanpy, duckdb; print('scdrs env ok')"
# wait (max 60 min) for the concurrently running LDSC job to finish writing comparator inputs
for i in $(seq 1 60); do [ -s data/processed/ldsc_in/BMI_Pulit2019.tsv ] && [ -s data/processed/ldsc_in/WHRadjBMI_Pulit2019.tsv ] && break; sleep 60; done
sleep 30; ls -la data/processed/ldsc_in/*Pulit2019.tsv 2>/dev/null
[ -s data/psd_meta_magma_scdrs.gs ] || PY=$ENVS/scdrs/bin/python bash scripts/13a_magma_gene_analysis.sh || { echo "MAGMA FAILED"; exit 1; }
$ENVS/scdrs/bin/python scripts/13b_scdrs_meta_magma.py || echo "SCDRS FAILED"
echo "=== magma+scdrs finished $(date) ==="
