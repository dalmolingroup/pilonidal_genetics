#!/bin/bash
# Submit the full pipeline as a dependency chain (run on the login node from the repo root).
# intel-256 jobs are moved to amd-512 (intel-256 nodes hung DuckDB jobs in Oct 2026).
set -euo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
sb() { sbatch --parsable "$@"; }
A=$(sb slurm/run_full_rerun_biallelic.sbatch)                                         # 01b-01e, 02a-02b
B=$(sb -p amd-512 -d afterok:$A slurm/run_magma.sbatch)                               # 03a
C=$(sb -p amd-512 -d afterok:$A slurm/run_shared_genetics_coloc_lava_gsem.sbatch)      # 04a-04d
D=$(sb -d afterok:$C slurm/run_eqtl_coloc.sbatch)                                      # 05a (needs 04a shared inputs)
F=$(sb -d afterok:$B slurm/run_magma_genesets.sbatch)                                  # 06a
I=$(sb -p intel-128 -c 4 --mem=32G -t 03:00:00 -J psd_tables -o slurm/psd_tables_%j.out -d afterok:$C:$D \
     --wrap "export NCPU=4; P=$ENVS/psd/bin/python; \$P scripts/06b_derived_summary_tables.py || echo DERIVED FAILED; \$P scripts/07a_locus_to_gene_table.py || echo L2G FAILED; \$P scripts/09a_lead_variant_supplementary_table.py || echo LEAD FAILED; \$P scripts/09b_fill_missing_rsids.py || echo RSID FAILED")
J=$(sb -d afterok:$I slurm/run_gwas_catalog.sbatch)                                    # 10a
echo "core=$A magma=$B shared=$C eqtl=$D genesets=$F tables=$I gwas_catalog=$J"
