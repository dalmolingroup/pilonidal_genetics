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
# Step 19a: supplementary lead-variant table (fixed and random effects, heterogeneity, per-cohort ORs, GRCh38 + GRCh37).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
zcat data/processed/meta/psd_meta_EUR.tsv.gz | head -1
NCPU=4 /mnt/beegfs/scratch/$USER/envs/psd/bin/python scripts/19a_lead_variant_supplementary_table.py || echo "LEAD TABLE FAILED"
