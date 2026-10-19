#!/usr/bin/env Rscript
# Step 14c: LAVA v0.1.5 local genetic correlation (Werme et al. 2022) across the 2,495 LD-independent genomic blocks,
# with the UK Biobank LAVA LD reference (v1.1, 100k unrelated Europeans; recommended over 1000G by the LAVA authors).
#  1) univariate local h2 test for every trait and block
#  2) bivariate local rg: PSD vs each trait, in blocks where both have univariate p < 1e-4
#  3) partial correlation PSD~HS | BMI in blocks with a nominally significant PSD-HS local rg where BMI also has local signal
# Sample overlap: standardized LDSC cross-trait intercepts (data/lava/sample.overlap.txt, built by the SLURM job).
# Outputs: results/lava_univariate.tsv, results/lava_bivariate_psd.tsv, results/lava_pcor_psd_hs_given_bmi.tsv
suppressMessages({library(LAVA); library(data.table); library(parallel)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
L <- file.path(base, "data/lava")
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
phenos <- c("PSD", "HS", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI")
input <- process.input(input.info.file = file.path(L, "input.info.txt"), sample.overlap.file = file.path(L, "sample.overlap.txt"),
                       ref.prefix = file.path(L, "ref", "lava-ukb-v1.1"), phenos = phenos)
loci <- read.loci(file.path(L, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
cat(nrow(loci), "loci;", ncores, "cores\n")
UNIV_P <- 1e-4

one_locus <- function(i) {
  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
  if (is.null(locus)) return(NULL)
  info <- data.table(locus = loci$LOC[i], chr = loci$CHR[i], start = loci$START[i], stop = loci$STOP[i], n_snps = locus$n.snps)
  u <- tryCatch(as.data.table(run.univ(locus)), error = function(e) NULL)
  if (is.null(u)) return(NULL)
  out <- list(univ = cbind(info, u), bivar = NULL, pcor = NULL)
  if (!("PSD" %in% u$phen) || u[phen == "PSD", p] >= UNIV_P) return(out)
  others <- u[phen != "PSD" & p < UNIV_P, phen]
  if (!length(others)) return(out)
  b <- tryCatch(as.data.table(run.bivar(locus, phenos = c("PSD", others), target = "PSD")), error = function(e) NULL)
  if (!is.null(b)) out$bivar <- cbind(info, b)
  if (!is.null(b) && "HS" %in% b$phen2 && b[phen2 == "HS", p] < 0.05 && "BMI" %in% others) {
    pc <- tryCatch(as.data.table(run.pcor(locus, target = c("PSD", "HS"), phenos = "BMI")), error = function(e) NULL)
    if (!is.null(pc)) out$pcor <- cbind(info, pc)
  }
  out
}
res <- mclapply(seq_len(nrow(loci)), one_locus, mc.cores = ncores, mc.preschedule = FALSE)
univ <- rbindlist(lapply(res, `[[`, "univ"), fill = TRUE)
bivar <- rbindlist(lapply(res, `[[`, "bivar"), fill = TRUE)
pcor <- rbindlist(lapply(res, `[[`, "pcor"), fill = TRUE)
if (nrow(bivar)) bivar[, p_bonferroni := pmin(1, p * .N)][, fdr := p.adjust(p, "BH")]
fwrite(univ, file.path(base, "results/lava_univariate.tsv"), sep = "\t")
fwrite(bivar, file.path(base, "results/lava_bivariate_psd.tsv"), sep = "\t")
fwrite(pcor, file.path(base, "results/lava_pcor_psd_hs_given_bmi.tsv"), sep = "\t")
cat("\nloci with PSD univariate p <", UNIV_P, ":", univ[phen == "PSD" & p < UNIV_P, .N], "\n")
cat("bivariate tests:", nrow(bivar), "; Bonferroni-significant:", bivar[p_bonferroni < 0.05, .N], "; FDR<0.05:", bivar[fdr < 0.05, .N], "\n")
if (nrow(bivar)) print(bivar[order(p)][1:min(.N, 40), .(locus, chr, start, stop, phen2, rho = round(rho, 3), rho.lower = round(rho.lower, 2), rho.upper = round(rho.upper, 2), p = signif(p, 3), fdr = signif(fdr, 3))])
if (nrow(pcor)) { cat("\npartial correlation PSD~HS | BMI:\n"); print(pcor) }
cat("LAVA_DONE\n")
