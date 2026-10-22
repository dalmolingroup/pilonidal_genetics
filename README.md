# Genetic architecture of pilonidal sinus disease (PSD)

A multi-biobank, multi-ancestry meta-analysis of PSD genome-wide association studies built entirely on open summary
statistics, and its overlap with the follicular occlusion disorders (hidradenitis suppurativa, acne), hair patterning
(balding) and adiposity.

## Main results

- **25 genome-wide significant loci outside the MHC, 19 novel** (6,771 EUR / 7,868 multi-ancestry cases). SNP heritability
  13–15% (liability scale); LDSC intercept 1.05.
- **PSD shares most of its common-variant liability with HS** (rg 0.76). The PSD–BMI correlation is explained by the
  HS-shared liability (Genomic SEM), and loci shared with BMI carry different causal variants (coloc).
- **A separate hair-patterning component**: PSD risk goes with balding genome-wide, in local blocks (LAVA) and at shared
  causal variants (*HDAC9*, *DKK1/MBL2*, *TCHH*, *SOX2*, chr20:22.1 Mb).
- **Candidate effector genes**: *SEMA4B* (skin eQTL), *WNT16* (fibroblast eQTL) and *WNT10A* (protective p.Phe228Ile).
- No robust cell-type or pathway signal once the MHC is removed.

All coordinates are **GRCh38** unless stated. MAGMA and LAVA run on GRCh37 references; every cross-dataset match is by
rsID, so builds are never mixed. The MHC (chr6:28–34 Mb) is excluded from locus counts, coloc, LAVA and gene sets.

## Data

| Role | Source |
|---|---|
| PSD GWAS | MVP multi-ancestry (GCST90480439) and EUR (GCST90478766), UK Biobank WGS EUR (GCST90473908), FinnGen R13 `L12_PILONIDALCYST`, Genes & Health SAS (GCST90727025) |
| Comparator traits | HS (GCST90726609; FinnGen R13 `L12_HIDRADENITISSUP`), acne (FinnGen R13 `L12_ACNE`; GCST90245818), male-pattern balding (GCST90043616), BMI and WHR adj. BMI (Pulit 2019) |
| eQTL | GTEx v8 skin, adipose and fibroblasts via the eQTL Catalogue |
| Single-cell atlases | CELLxGENE skin atlases, GSE274880 (HS skin), GSE287304 (hair-follicle mesenchyme) |
| References | 1000G EUR (MAGMA, LDSC), UKB LAVA reference v1.1, MSigDB 2025.1, Open Targets credible sets |

Raw data and large intermediates live on the cluster under `data/raw/` and `data/processed/` (git-ignored).

## Pipeline

Heavy steps run on the NPAD SLURM cluster (`sc2.npad.ufrn.br`, repository at `~/scratch/pilonidal_genetics`); the batch
scripts are in `slurm/`. Steps marked *local* are lightweight.

| Step | Script | What it does | Main outputs (`results/`) |
|---|---|---|---|
| 01a | `01a_harmonize_psd_cohorts.sh` | Harmonise the 5 PSD cohorts to GRCh38 | `data/processed/harmonized/` |
| 01b | `01b_meta_analyze_psd_gwas.py` | IVW fixed-effect meta-analyses META_EUR and META_ALL, cohort QC | `psd_meta_{EUR,ALL}_loci.tsv`, `psd_meta_cohort_qc.tsv`, `psd_published_loci_replication.tsv` |
| 01c | `01c_finemap_psd_loci.py` | ABF and SuSiE fine-mapping, VEP and GTEx annotation | `psd_finemap_*.tsv` |
| 01d | `01d_classify_meta_loci.py` *(local)* | Merge loci, known vs novel, intercept correction, leave-one-cohort-out | `psd_meta_locus_classification.tsv`, `psd_meta_loci_intercept_corrected.tsv`, `psd_meta_loci_nonfinngen_replication.tsv` |
| 01e | `01e_reannotate_vep_by_rsid.py` *(local)* | Re-annotate credible-set variants with VEP by rsID | `psd_finemap_vep_by_rsid.tsv` |
| 02a–b | `02a_prepare_ldsc_inputs.py`, `02b_run_ldsc_rg.sh` | LDSC h2, intercepts and genetic correlations | `ldsc_v2_*.tsv`, `ldsc_logs_v2/` |
| 03a–c | `03a_magma_gene_analysis.sh`, `03b_scdrs_meta_magma.py`, `03c_magma_genesets_exclude_mhc.sh` | MAGMA gene scores; scDRS cell-type mapping with and without the MHC | `magma_genes_*.tsv`, `scdrs_magma_*.tsv` |
| 04a–d | `04a_prepare_shared_genetics_inputs.py`, `04b_coloc_psd_vs_traits.R`, `04c_lava_local_rg.R`, `04d_genomic_sem.R` | Shared genetics with HS, acne, balding, BMI, WHR: coloc, LAVA local rg, Genomic SEM | `coloc_psd_vs_traits.tsv`, `lava_bivariate_psd_noMHC.tsv`, `genomicsem_*.tsv`, `hs_gws_loci.tsv` |
| 05a | `05a_eqtl_coloc_gtex.R` | eQTL colocalisation in GTEx skin, adipose and fibroblasts | `eqtl_coloc_psd_gtex.tsv` |
| 06a | `06a_magma_geneset_tissue.sh` | MHC-free MAGMA gene-set (GO BP, Reactome) and GTEx tissue analyses | `magma_genesets_*.tsv`, `magma_gtex_tissues_nomhc.tsv` |
| 07a | `07a_locus_to_gene_table.py` *(local)* | Locus-to-gene evidence scoring for the 25 loci | `psd_locus_to_gene.tsv` |
| 08a–b | `08a_network_propagation.py`, `08b_netprop_chromosome_split.py` | Network propagation (Visonà et al. 2024) and its chromosome-split validation | `netprop_*.tsv` |
| 09a | `09a_lead_variant_supplementary_table.py` | Lead variants: fixed and random effects, heterogeneity, per-cohort ORs, GRCh38 and GRCh37 positions | `psd_lead_variants_supplementary.tsv` |
| 09b | `09b_fill_missing_rsids.py` *(local)* | Fill indel rsIDs via Ensembl VEP; rerun after regenerating steps 01–05 | updates the tables above |

`scripts/lib/` holds two helper modules, `scdrs_atlases.py` and `gse274880_hs_atlas.py`, that step 03b imports for atlas
download and annotation. They come from the early exploratory pipeline, which was removed.

Notes:
- Network propagation did not pass its pre-set validation (it did not beat degree-matched permutations), so it is not
  used in the results.
- The 4 loci with significant between-cohort heterogeneity (*DKK1/MBL2*, *TCHH*, *KDM1A*, *HDAC9*) are flagged in
  `psd_lead_variants_supplementary.tsv`, with random-effects estimates as a sensitivity analysis.

## Repository layout

- `scripts/`: the pipeline above, numbered in run order; `scripts/lib/` holds shared helpers.
- `slurm/`: batch scripts and the input download script (`cluster_download_inputs.sh`).
- `data/`: small inputs kept in git (`psd_known_loci_prior_art.tsv`, used by step 01d, and the MAGMA scDRS gene sets `psd_meta_magma_scdrs.gs`); everything else is git-ignored.
- `results/`: result tables; `results/cluster_pull/` holds the raw pulls and job logs from the cluster runs.
- `figures/`: meta-analysis Manhattan and QQ plot.
- `presentation/`: Quarto talk (`psd_genetics_talk.qmd`), its figures (`make_figures.py`) and references.
- `resources/`: papers consulted (not tracked).
