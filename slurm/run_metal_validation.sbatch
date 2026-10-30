#!/bin/bash
#SBATCH --job-name=psd_metal
#SBATCH --partition=amd-512
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/psd_metal_%j.out
#SBATCH --error=slurm/psd_metal_%j.err
# Step 01f: reproduce META_ALL and META_EUR with METAL and compare with step 01b.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
LOCAL=/scratch/local/$USER/job_$SLURM_JOB_ID; mkdir -p $LOCAL/duckdb; trap 'rm -rf $LOCAL' EXIT
cp -a $ENVS/psd $LOCAL/psd
# METAL 2011-03-25 Linux binary
[ -x $ENVS/metal/generic-metal/metal ] || { mkdir -p $ENVS/metal && curl -fsSL -o $ENVS/metal/metal.tgz \
  https://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz && tar xzf $ENVS/metal/metal.tgz -C $ENVS/metal; }
export METAL=$ENVS/metal/generic-metal/metal DUCKDB_MEM=80GB DUCKDB_THREADS=16 DUCKDB_TMP=$LOCAL/duckdb
$LOCAL/psd/bin/python scripts/01f_metal_validation.py
