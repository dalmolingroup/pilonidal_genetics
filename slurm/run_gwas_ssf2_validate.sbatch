#!/bin/bash
#SBATCH --job-name=psd_ssf2
#SBATCH --partition=intel-128
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/psd_ssf2_%j.out
#SBATCH --error=slurm/psd_ssf2_%j.err
# Step 10a check: re-validate the GWAS-SSF files with gwas-sumstats-tools 2.0.2 (needs Python >= 3.13).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
E=/mnt/beegfs/scratch/$USER/envs/gwas_ssf2
module load softwares/python/3.13.7-gnu8
[ -x $E/bin/gwas-ssf ] || { python3 -m venv $E && $E/bin/pip install -q "gwas-sumstats-tools==2.0.2"; }
$E/bin/pip show gwas-sumstats-tools | head -2
(cd data/processed/gwas_catalog && md5sum -c md5sums.txt)
for f in data/processed/gwas_catalog/psd_meta_*_GRCh38.tsv.gz; do echo "== validate $f"; $E/bin/gwas-ssf validate "$f" 2>&1 | tail -15; echo "exit ${PIPESTATUS[0]}"; done
