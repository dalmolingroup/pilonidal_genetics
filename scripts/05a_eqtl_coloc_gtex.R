#!/usr/bin/env Rscript
# Step 05a: Bayesian colocalization (coloc.abf; priors as step 04b) of the non-MHC PSD loci with GTEx v8 cis-eQTLs, using the
# eQTL Catalogue uniformly processed summary statistics: skin sun-exposed (QTD000316), skin not sun-exposed (QTD000311),
# subcutaneous adipose (QTD000116), visceral adipose (QTD000121) and cultured fibroblasts (QTD000216).
# Region = PSD lead +/-500 kb (GRCh38); variants matched by position + alleles; a gene is tested when its eQTL p < 1e-4 in the window.
# Inputs : data/processed/eqtl/<dataset>/chr<chr>_<pos>.tsv.gz (tabix extracts from slurm/run_eqtl_coloc.sbatch),
#          data/processed/eqtl/ensembl105_genes.tsv, data/processed/shared/PSD_META_{EUR,ALL}.tsv.gz (step 04a)
# Output : results/eqtl_coloc_psd_gtex.tsv
suppressMessages({library(data.table); library(coloc)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
eq <- file.path(base, "data/processed/eqtl")
datasets <- data.table(dataset = c("QTD000316", "QTD000311", "QTD000116", "QTD000121", "QTD000216"),
                       tissue = c("Skin_sun_exposed", "Skin_not_sun_exposed", "Adipose_subcutaneous", "Adipose_visceral", "Fibroblast_cultured"))
cols <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type",
          "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
genes <- fread(file.path(eq, "ensembl105_genes.tsv"), header = FALSE, col.names = c("gene_id", "gene_name")); setkey(genes, gene_id)
loci <- fread(file.path(base, "results/psd_meta_locus_classification.tsv"))[as.character(is_mhc) %in% c("False", "FALSE")]
loci[, label := ifelse(is.na(nearest_gene) | nearest_gene == "", paste0("chr", chr, ":", pos), nearest_gene)]
loci[, tag := paste0("chr", chr, "_", pos)]
psd_sets <- list(PSD_META_EUR = 6771 / 1381000, PSD_META_ALL = 7868 / 1606360)
cat(nrow(loci), "non-MHC PSD loci x", nrow(datasets), "GTEx tissues\n")

out <- list()
for (pset in names(psd_sets)) {
  psd <- fread(file.path(base, "data/processed/shared", paste0(pset, ".tsv.gz")), showProgress = FALSE); psd[, CHR := as.character(CHR)]
  for (i in seq_len(nrow(loci))) {
    L <- loci[i]
    w <- psd[CHR == as.character(L$chr) & POS >= L$pos - 5e5 & POS <= L$pos + 5e5]
    if (nrow(w) < 50) next
    for (j in seq_len(nrow(datasets))) {
      f <- file.path(eq, datasets$dataset[j], paste0(L$tag, ".tsv.gz"))
      if (!file.exists(f) || file.size(f) < 100) { cat("  missing", f, "\n"); next }
      e <- fread(f, header = FALSE, col.names = cols, showProgress = FALSE)[type == "SNP" & se > 0 & maf > 0 & is.finite(beta)]
      gtest <- e[, .(minp = min(pvalue)), by = molecular_trait_id][minp < 1e-4, molecular_trait_id]
      for (g in gtest) {
        m <- merge(w[, .(POS, SNP, A1, A2, b1 = BETA, se1 = SE, p1 = P, eaf1 = EAF, N1 = N)],
                   e[molecular_trait_id == g, .(POS = position, ref, alt, b2 = beta, se2 = se, p2 = pvalue, maf2 = maf, an)], by = "POS")
        m <- m[(A1 == alt & A2 == ref) | (A1 == ref & A2 == alt)]
        m[A1 == ref & A2 == alt, b2 := -b2]                                  # eQTL beta is for alt; align to the PSD effect allele A1
        m <- m[!duplicated(SNP)]
        if (nrow(m) < 50) next
        d1 <- list(beta = m$b1, varbeta = m$se1^2, snp = m$SNP, type = "cc", s = psd_sets[[pset]], N = median(m$N1), MAF = pmin(m$eaf1, 1 - m$eaf1))
        d2 <- list(beta = m$b2, varbeta = m$se2^2, snp = m$SNP, type = "quant", N = median(m$an) / 2, MAF = m$maf2)
        res <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)), error = function(e) NULL)
        res2 <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)), error = function(e) NULL)
        if (is.null(res)) next
        s <- as.list(res$summary); top <- res$results[which.max(res$results$SNP.PP.H4), ]
        jt <- m[SNP == top$snp][1]
        out[[length(out) + 1]] <- data.table(psd_set = pset, label = L$label, chr = L$chr, pos = L$pos, psd_lead = L$rsid, psd_status = L$status,
          dataset = datasets$dataset[j], tissue = datasets$tissue[j], gene_id = g, gene_name = genes[J(g), gene_name],
          nsnps = s$nsnps, PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf, PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf,
          PP.H4_p12_1e5 = if (!is.null(res2)) res2$summary[["PP.H4.abf"]] else NA, min_p_psd = min(m$p1), min_p_eqtl = min(m$p2),
          h4_top_snp = top$snp, h4_top_snp_pp = top$SNP.PP.H4, h4_top_beta_psd = jt$b1, h4_top_beta_eqtl = jt$b2,
          risk_allele_expression = ifelse(sign(jt$b1) * jt$b2 > 0, "up", "down"))
      }
    }
    cat(pset, L$label, "done\n")
  }
}
res <- rbindlist(out)
fwrite(res, file.path(base, "results/eqtl_coloc_psd_gtex.tsv"), sep = "\t")
cat("\n", nrow(res), "locus x gene x tissue tests\n")
cat("\n=== PP.H4 >= 0.8 (shared causal variant), PSD_META_EUR ===\n")
print(res[psd_set == "PSD_META_EUR" & PP.H4 >= 0.8, .(label, tissue, gene_name, nsnps, PP.H4 = round(PP.H4, 3), PP.H3 = round(PP.H3, 3),
          h4_top_snp, h4_top_snp_pp = round(h4_top_snp_pp, 2), risk_allele_expression)][order(label, -PP.H4)])
cat("\n=== 0.5 <= PP.H4 < 0.8 (suggestive), PSD_META_EUR ===\n")
print(res[psd_set == "PSD_META_EUR" & PP.H4 >= 0.5 & PP.H4 < 0.8, .(label, tissue, gene_name, PP.H4 = round(PP.H4, 3), PP.H3 = round(PP.H3, 3))][order(label, -PP.H4)])
cat("EQTL_COLOC_DONE\n")
