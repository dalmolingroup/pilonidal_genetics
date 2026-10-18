# Pilonidal Sinus Disease (PSD) Multi-Omics & Single-Cell Pathogenesis

Computational dissection of Pilonidal Sinus Disease (ICD-10 `L05` / PheCode `686.3`) genetic risk, cell-type specificity, and shared/divergent architecture with follicular occlusion and ectodermal/mesenchymal patterning traits.

## Directions
1. **Direction 1 — Single-Cell & Regulatory Dissection of PSD Risk Loci (`scDRS` + Cell-Type Specificity + eQTL Colocalization)**
   - Harmonize open GWAS/exome summary statistics (`GCST90480439`, `GCST90478766`, `GCST90473908`, `GCST90727025`, `GCST90080332`, `GCST90084318`) and Open Targets L2G/genetic associations (`HP_0010769`, `HP_0010771`, `MONDO_0008249`).
   - Map polygenic and locus-prioritized gene sets onto human skin, hair follicle, and hidradenitis suppurativa single-cell atlases via `scDRS` and cell-type specificity scoring.
   - Colocalize non-coding risk loci (*HDAC9–TWIST1*, *DKK1–MBL2*, *TBX15–WARS2*, *PRKAR1B*, *NEK9*, *CHRM3*, *RNFT2*) with GTEx skin/adipose eQTLs/sQTLs and single-cell regulatory compartments.
2. **Direction 2 — Multi-Ancestry Meta-GWAS & Fine-Mapping (EUR + AFR + AMR + SAS)**
3. **Direction 3 — Genetic & Transcriptomic Architecture of the Follicular Occlusion Tetrad (PSD vs. Hidradenitis Suppurativa vs. Hair/Adipose Patterning)**

## Repository Structure
- `data/`: Metadata tables, gene sets, and pointers/scripts for downloading large GWAS summary statistics and `.h5ad` atlases (heavy raw files live on SLURM `scratch/` or ignored locally).
- `scripts/`: Modular Python/R scripts for GWAS harmonization, gene-level scoring (`MAGMA`/locus mapping), `scDRS`, and eQTL colocalization.
- `slurm/`: Batch scripts configured for `sc2.npad.ufrn.br`.
- `results/`: Summary tables, enrichment statistics, and colocalization outputs.
- `figures/`: Publication-grade visualizations.
