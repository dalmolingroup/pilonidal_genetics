#!/bin/bash
#SBATCH --job-name=psd_eqtl_coloc
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/psd_eqtl_coloc_%j.out
#SBATCH --error=slurm/psd_eqtl_coloc_%j.err
# Step 15a: GTEx v8 eQTL colocalization of the non-MHC PSD loci.
# (1) remote tabix extracts (lead +/-500 kb) of eQTL Catalogue summary statistics for 5 GTEx tissues, one process per dataset;
# (2) Ensembl 105 gene id -> name table; (3) scripts/15a_eqtl_coloc_gtex.R
set -uo pipefail
cd /mnt/beegfs/scratch/$USER/pilonidal_genetics
ENVS=/mnt/beegfs/scratch/$USER/envs
export PATH=$ENVS/rgen/bin:$PATH
EQ=data/processed/eqtl; mkdir -p $EQ
TABIX=$(command -v tabix || true)
if [ -z "$TABIX" ]; then
  CONDA=/opt/npad/shared/conda/forge2023/bin
  [ -x $ENVS/htslib/bin/tabix ] || $CONDA/mamba create -y -q -p $ENVS/htslib -c conda-forge -c bioconda htslib > /dev/null 2>&1 \
    || $CONDA/conda create -y -q -p $ENVS/htslib -c conda-forge -c bioconda htslib > /dev/null 2>&1
  TABIX=$ENVS/htslib/bin/tabix
fi
$TABIX --version | head -1 || { echo "NO TABIX"; exit 1; }
export TABIX
# non-MHC PSD loci (is_mhc is the last column): chr pos window_start window_end
awk -F'\t' 'NR>1 && $NF=="False"{s=$2-500000; if (s<1) s=1; print $1, $2, s, $2+500000}' results/psd_meta_locus_classification.tsv > $EQ/loci.txt
echo "$(wc -l < $EQ/loci.txt) loci"
fetch() {  # $1 dataset; the remote .tbi is cached in the dataset directory, so one process per dataset
  ds=$1; d=data/processed/eqtl/$ds; mkdir -p $d; cd $d || exit 1
  url=https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/$ds/$ds.all.tsv.gz
  while read chr pos s e; do
    out=chr${chr}_${pos}.tsv.gz
    [ -s $out ] && continue
    for try in 1 2 3; do $TABIX $url $chr:$s-$e | gzip > $out.tmp && mv $out.tmp $out && break; sleep 15; done
  done < ../loci.txt
  echo "$ds: $(ls chr*.tsv.gz 2>/dev/null | wc -l) loci extracted"
}
export -f fetch
printf "QTD000316\nQTD000311\nQTD000116\nQTD000121\nQTD000216\n" | xargs -P 5 -I{} bash -c 'fetch {}'
[ -s $EQ/ensembl105_genes.tsv ] || curl -sSL --retry 5 https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz | zcat \
  | awk -F'\t' '$3=="gene"{id=$9; sub(/.*gene_id "/,"",id); sub(/".*/,"",id); nm=id; if ($9 ~ /gene_name "/) {nm=$9; sub(/.*gene_name "/,"",nm); sub(/".*/,"",nm)}; print id"\t"nm}' > $EQ/ensembl105_genes.tsv
echo "$(wc -l < $EQ/ensembl105_genes.tsv) Ensembl genes"
Rscript scripts/15a_eqtl_coloc_gtex.R || echo "EQTL COLOC FAILED"
echo "=== eqtl coloc finished $(date) ==="
