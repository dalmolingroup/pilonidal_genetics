#!/bin/bash
#SBATCH --job-name=psd_leads
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm/psd_leads_%j.out
#SBATCH --error=slurm/psd_leads_%j.err
# Step 09a: supplementary lead-variant table (fixed and random effects, heterogeneity, per-cohort ORs, GRCh38 + GRCh37).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
module load softwares/R/4.5.0-gnu11.2 2>/dev/null || true
zcat data/processed/meta/psd_meta_EUR.tsv.gz | head -1
NCPU=4 METAFOR_RSCRIPT=${METAFOR_RSCRIPT:-Rscript} /mnt/beegfs/scratch/$USER/envs/psd/bin/python scripts/09a_lead_variant_supplementary_table.py || echo "LEAD TABLE FAILED"
