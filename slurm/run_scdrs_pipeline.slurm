#!/bin/bash
#SBATCH --job-name=psd_scdrs
#SBATCH --output=slurm/psd_scdrs_%j.out
#SBATCH --error=slurm/psd_scdrs_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# Assumes repository is cloned into /scratch/$USER/pilonidal_genetics on sc2.npad.ufrn.br
WORKDIR="${SLURM_SUBMIT_DIR:-/scratch/$USER/pilonidal_genetics}"
cd "$WORKDIR"

if [ ! -d ".venv" ]; then
    uv venv .venv
    uv pip install --python .venv/bin/python scdrs scanpy anndata pandas scipy statsmodels matplotlib seaborn
fi

source .venv/bin/activate

python scripts/01_fetch_genetics_and_exome.py
python scripts/02_stream_gwas_and_build_genesets.py
python scripts/03_run_scdrs_and_specificity.py
python scripts/04_gtex_eqtl_and_hair_follicle_mesenchyme.py
