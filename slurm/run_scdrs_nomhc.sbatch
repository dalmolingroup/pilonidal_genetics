#!/bin/bash
#SBATCH --job-name=psd_scdrs_nomhc
#SBATCH --partition=intel-256
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --array=0-6
#SBATCH --output=slurm/psd_scdrs_nomhc_%A_%a.out
#SBATCH --error=slurm/psd_scdrs_nomhc_%A_%a.err
# scDRS with MHC-excluded MAGMA gene sets (step 03c), one atlas per array task. Merge afterwards with SCDRS_MERGE=1 (see bottom).
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2   # must match the values script 09 sets at import
module load softwares/python/3.13.7-gnu8   # venv python links against this module's libpython3.13.so
export PATH=$ENVS/scdrs/bin:$PATH
ATLASES=(tabula_sapiens_skin skin_12datasets_integrated skin_keratinocyte_subclusters skin_fibroblast_subclusters
         scalp_trunk_skin_68k GSE274880_HS GSE287304_hair_follicle_mesenchyme)
GS=data/psd_meta_magma_scdrs_noMHC.gs
# task 0 builds the gene sets; the others wait for it
if [ "$SLURM_ARRAY_TASK_ID" = 0 ]; then
  [ -s $GS ] || PY=$ENVS/scdrs/bin/python bash scripts/03c_magma_genesets_exclude_mhc.sh || { echo "GENESETS FAILED"; exit 1; }
else
  for i in $(seq 1 30); do [ -s $GS ] && break; sleep 20; done
fi
sleep 5
SCDRS_GS=$GS SCDRS_TAG=scdrs_magma_nomhc SCDRS_ATLAS=${ATLASES[$SLURM_ARRAY_TASK_ID]} $ENVS/scdrs/bin/python scripts/03b_scdrs_meta_magma.py || echo "SCDRS FAILED"
echo "=== ${ATLASES[$SLURM_ARRAY_TASK_ID]} finished $(date) ==="
# merge (dependent job): SCDRS_TAG=scdrs_magma_nomhc SCDRS_MERGE=1 python scripts/03b_scdrs_meta_magma.py
