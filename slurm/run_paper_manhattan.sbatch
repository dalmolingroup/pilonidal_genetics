#!/bin/bash
#SBATCH --job-name=psd_manhattan
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=slurm/psd_manhattan_%j.out
#SBATCH --error=slurm/psd_manhattan_%j.err
# Paper Figure 1 Miami plot and supplementary QQ plots from the full meta-analysis summary statistics.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
NCPU=4 /mnt/beegfs/scratch/$USER/envs/psd/bin/python paper/cluster_manhattan.py || echo "MANHATTAN FAILED"
