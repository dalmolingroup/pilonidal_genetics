#!/bin/bash
#SBATCH --job-name=psd_netprop
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/psd_netprop_%j.out
#SBATCH --error=slurm/psd_netprop_%j.err
# Step 18a: network propagation (personalized PageRank, Visona et al. 2024) of MHC-free MAGMA gene scores over 5 networks.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
module load softwares/python/3.13.7-gnu8   # scdrs venv python links against this module's libpython
NCPU=5 $ENVS/scdrs/bin/python scripts/18a_network_propagation.py || echo "NETPROP FAILED"
echo "=== network propagation finished $(date) ==="
