#!/bin/bash
#SBATCH --job-name=psd_netprop_val
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/psd_netprop_val_%j.out
#SBATCH --error=slurm/psd_netprop_val_%j.err
# Step 18b: chromosome-split validation of the network propagation (go/no-go vs degree baseline, 1,000 degree-matched permutations).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
module load softwares/python/3.13.7-gnu8   # scdrs venv python links against this module's libpython
NCPU=20 $ENVS/scdrs/bin/python scripts/18b_netprop_chromosome_split.py || echo "NETPROP VALIDATION FAILED"
echo "=== network propagation validation finished $(date) ==="
