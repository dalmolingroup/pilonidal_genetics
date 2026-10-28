#!/bin/bash
#SBATCH --job-name=psd_full_rerun
#SBATCH --partition=amd-512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/psd_full_rerun_%j.out
#SBATCH --error=slurm/psd_full_rerun_%j.err
# Re-run after excluding multi-allelic sites (01b) and making UKB LD the primary SuSiE reference with coloc::finemap.abf (01c):
# meta -> fine-mapping -> locus classification -> VEP by rsID -> LDSC. Downstream steps are chained by slurm/submit_full_rerun.sh.
# amd-512: intel-256 nodes hung DuckDB queries for hours (Oct 2026). DuckDB spills to node-local disk, not beegfs.
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
BASE=$PWD; ENVS=/mnt/beegfs/scratch/$USER/envs
# Run the Python env from node-local disk: processes stalled in the kernel waiting on beegfs-backed pages of the env's libraries.
LOCAL=/scratch/local/$USER/job_$SLURM_JOB_ID; mkdir -p $LOCAL/duckdb; trap 'rm -rf $LOCAL' EXIT
cp -a $ENVS/psd $LOCAL/psd
export DUCKDB_MEM=150GB DUCKDB_THREADS=${SLURM_CPUS_PER_TASK:-32} DUCKDB_TMP=$LOCAL/duckdb OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export PLINK2=$BASE/bin/plink2 LDSC_ENV=$ENVS/ldsc COLOC_RSCRIPT=$ENVS/rgen/bin/Rscript
export PATH=$ENVS/rsusie/bin:$PATH
PY=$LOCAL/psd/bin/python
$PY -c "import duckdb, sys; print(sys.prefix, duckdb.__file__)"
Rscript -e 'cat("susieR", requireNamespace("susieR",quietly=TRUE), "\n")'
$COLOC_RSCRIPT -e 'cat("coloc", as.character(packageVersion("coloc")), "\n")'

# back up the current results, then clear every cache that depends on the PSD meta-analysis
BK=results_backup_pre_biallelic
[ -d $BK ] || { cp -a results $BK && mkdir -p $BK/meta && cp -a data/processed/meta/psd_meta_*.tsv.gz $BK/meta/; }
rm -f results/psd_meta_ALL_loci.tsv results/psd_meta_EUR_loci.tsv results/psd_finemap_locus_summary.tsv
rm -rf data/processed/ldsc_in data/processed/ldsc_munged results/ldsc_logs results/ldsc_overlap
rm -f data/magma/in/PSD_* data/magma/out/PSD_* data/magma/out/magma_zscores* data/magma/out_nomhc/PSD_* data/psd_meta_magma_scdrs*.gs
rm -f data/processed/shared/PSD_* data/lava/sumstats/PSD_* data/netprop/PSD_META_*.npz
rm -f results/coloc_psd_vs_traits.tsv results/genomicsem_models.tsv results/lava_bivariate_psd.tsv

$PY scripts/01b_meta_analyze_psd_gwas.py || { echo "META FAILED"; exit 1; }
$PY scripts/01c_finemap_psd_loci.py || echo "FINEMAP FAILED"
$PY scripts/01d_classify_meta_loci.py || echo "CLASSIFY FAILED"
$PY scripts/01e_reannotate_vep_by_rsid.py || echo "VEP FAILED"
$PY scripts/02a_prepare_ldsc_inputs.py && bash scripts/02b_run_ldsc_rg.sh || echo "LDSC FAILED"
echo "=== full rerun core finished $(date) ==="
