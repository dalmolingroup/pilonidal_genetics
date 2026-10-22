#!/bin/bash
#SBATCH --job-name=psd_susie
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/psd_susie_%j.out
#SBATCH --error=slurm/psd_susie_%j.err
# Re-run fine-mapping (01c) with SuSiE available: susieR from conda-forge (module R lacks a matching libstdc++).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD; ENVS=/mnt/beegfs/scratch/$USER/envs
export PLINK2=$BASE/bin/plink2 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
if [ ! -x $ENVS/rsusie/bin/Rscript ]; then
  echo "[env] creating conda env with r-susier"; ~/miniconda3/bin/conda create -y -q -p $ENVS/rsusie -c conda-forge r-base=4.3 r-susier r-data.table > slurm/rsusie_env_create.log 2>&1 || echo "[env] r-susier env creation FAILED"
fi
export PATH=$ENVS/rsusie/bin:$PATH
Rscript -e 'cat("susieR", requireNamespace("susieR",quietly=TRUE), as.character(packageVersion("susieR")), "\n")'
if [ ! -s results/psd_meta_ALL_loci.tsv ]; then echo "meta loci missing; main job did not finish"; exit 1; fi
mkdir -p results/finemap_abf_only && cp -f results/psd_finemap_*.tsv results/finemap_abf_only/ 2>/dev/null || true
rm -f results/psd_finemap_locus_summary.tsv
$ENVS/psd/bin/python scripts/01c_finemap_psd_loci.py
echo "=== susie rerun finished $(date) ==="
