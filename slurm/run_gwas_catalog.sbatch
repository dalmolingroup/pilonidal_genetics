#!/bin/bash
#SBATCH --job-name=psd_gwascat
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/psd_gwascat_%j.out
#SBATCH --error=slurm/psd_gwascat_%j.err
# Step 10a: GWAS-SSF files for GWAS Catalog deposition, then the official validator (gwas-sumstats-tools, own venv).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
NCPU=4 /mnt/beegfs/scratch/$USER/envs/psd/bin/python scripts/10a_format_gwas_catalog.py || { echo "FORMAT FAILED"; exit 1; }
for f in data/processed/gwas_catalog/psd_meta_*_GRCh38.tsv.gz; do
    echo "== validate $f"
    /mnt/beegfs/scratch/$USER/envs/gwas_ssf/bin/gwas-ssf validate "$f" || echo "VALIDATION FAILED: $f"
done
