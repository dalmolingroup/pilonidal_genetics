#!/usr/bin/env Rscript
# Step 04b: Bayesian colocalization (coloc.abf; Giambartolomei 2014, single causal variant per trait per region)
# of PSD with HS, acne, balding, BMI and WHRadjBMI at (i) all 26 PSD genome-wide significant loci and (ii) all HS
# genome-wide significant loci. Region = +/-500 kb around the lead (GRCh38, PSD coordinates); variants matched by rsID
# with allele alignment. Priors p1 = p2 = 1e-4; p12 = 5e-6 (Wallace 2020) with p12 = 1e-5 as sensitivity.
# Outputs: results/coloc_psd_vs_traits.tsv
suppressMessages({library(data.table); library(coloc)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
sh <- file.path(base, "data/processed/shared")
rd <- function(t) { d <- fread(file.path(sh, paste0(t, ".tsv.gz")), showProgress = FALSE); d[, CHR := as.character(CHR)]; setkey(d, SNP); d }
traits <- list(
  HS        = list(type = "cc",    s = 3941 / 1439544),
  HS_FG     = list(type = "cc",    s = 1453 / 476349),
  ACNE_FG   = list(type = "cc",    s = 4812 / 479708),
  ACNE_TL   = list(type = "cc",    s = 34422 / 399413),
  BALDP1    = list(type = "cc",    s = 66172 / 207036),
  BMI       = list(type = "quant", s = NA),
  WHRadjBMI = list(type = "quant", s = NA))
psd_sets <- list(PSD_META_EUR = 6771 / 1381000, PSD_META_ALL = 7868 / 1606360)

loci_psd <- fread(file.path(base, "results/psd_meta_locus_classification.tsv"))
regions <- loci_psd[, .(region_source = "PSD", chr = as.character(chr), pos = as.integer(pos), lead = rsid,
                        label = ifelse(is.na(nearest_gene) | nearest_gene == "", paste0("chr", chr, ":", pos), nearest_gene), psd_status = status)]
hs_loci <- fread(file.path(base, "results/hs_gws_loci.tsv"))
hs_reg <- hs_loci[, .(region_source = "HS", chr = as.character(CHR), pos = as.integer(POS_b38), lead = SNP, label = paste0("HS_chr", CHR, ":", POS_b38), psd_status = NA_character_)]
# drop HS regions already covered by a PSD region
hs_reg <- hs_reg[!sapply(seq_len(nrow(hs_reg)), function(i) any(regions$chr == hs_reg$chr[i] & abs(regions$pos - hs_reg$pos[i]) < 500000))]
regions <- rbind(regions, hs_reg)
cat(nrow(regions), "regions (", sum(regions$region_source == "PSD"), "PSD,", sum(regions$region_source == "HS"), "HS-only )\n")

tr <- lapply(names(traits), rd); names(tr) <- names(traits)
out <- list()
for (pset in names(psd_sets)) {
  psd <- rd(pset)
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    w <- psd[CHR == r$chr & POS >= r$pos - 5e5 & POS <= r$pos + 5e5]
    if (nrow(w) < 50) next
    for (t in names(traits)) {
      o <- tr[[t]][J(w$SNP), nomatch = 0L]
      m <- merge(w[, .(SNP, A1, A2, b1 = BETA, se1 = SE, p1 = P, eaf1 = EAF, N1 = N)], o[, .(SNP, oA1 = A1, oA2 = A2, b2 = BETA, se2 = SE, p2 = P, eaf2 = EAF, N2 = N)], by = "SNP")
      m <- m[(A1 == oA1 & A2 == oA2) | (A1 == oA2 & A2 == oA1)]
      m[A1 == oA2 & A2 == oA1, `:=`(b2 = -b2, eaf2 = 1 - eaf2)]
      if (nrow(m) < 50) next
      d1 <- list(beta = m$b1, varbeta = m$se1^2, snp = m$SNP, type = "cc", s = psd_sets[[pset]], N = median(m$N1), MAF = pmin(m$eaf1, 1 - m$eaf1))
      d2 <- list(beta = m$b2, varbeta = m$se2^2, snp = m$SNP, type = traits[[t]]$type, N = median(m$N2), MAF = pmin(m$eaf2, 1 - m$eaf2))
      if (traits[[t]]$type == "cc") d2$s <- traits[[t]]$s else d2$sdY <- 1
      res <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)), error = function(e) NULL)
      res2 <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)), error = function(e) NULL)
      if (is.null(res)) next
      s <- as.list(res$summary); top <- res$results[which.max(res$results$SNP.PP.H4), ]
      jt <- m[SNP == top$snp]
      out[[length(out) + 1]] <- data.table(psd_set = pset, region_source = r$region_source, label = r$label, chr = r$chr, pos = r$pos, psd_lead = r$lead,
        psd_status = r$psd_status, trait = t, nsnps = s$nsnps, PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf, PP.H3 = s$PP.H3.abf,
        PP.H4 = s$PP.H4.abf, PP.H4_p12_1e5 = if (!is.null(res2)) res2$summary[["PP.H4.abf"]] else NA,
        min_p_psd = min(m$p1), min_p_trait = min(m$p2), h4_top_snp = top$snp, h4_top_snp_pp = top$SNP.PP.H4,
        h4_top_beta_psd = jt$b1, h4_top_beta_trait = jt$b2, same_direction = sign(jt$b1) == sign(jt$b2))
    }
    cat(pset, i, r$label, "done\n")
  }
}
res <- rbindlist(out)
fwrite(res, file.path(base, "results/coloc_psd_vs_traits.tsv"), sep = "\t")
cat("\n=== PP.H4 >= 0.8 (strong colocalization), PSD_META_EUR ===\n")
print(res[psd_set == "PSD_META_EUR" & PP.H4 >= 0.8, .(label, trait, nsnps, PP.H4 = round(PP.H4, 3), PP.H3 = round(PP.H3, 3), h4_top_snp, h4_top_snp_pp = round(h4_top_snp_pp, 2), same_direction)][order(label)])
cat("\n=== PP.H3 >= 0.8 (distinct causal variants), PSD_META_EUR ===\n")
print(res[psd_set == "PSD_META_EUR" & PP.H3 >= 0.8, .(label, trait, PP.H3 = round(PP.H3, 3), min_p_psd, min_p_trait)][order(label)])
cat("COLOC_DONE\n")
