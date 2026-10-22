#!/usr/bin/env Rscript
# Step 04d: Genomic SEM (Grotzinger et al. 2019). Multivariable LDSC on HapMap3-munged summary statistics, then:
#  M1  common "follicular occlusion" factor  F =~ PSD + HS + ACNE_FG + ACNE_TL   (how much PSD genetic variance is shared vs PSD-specific)
#  M1b same factor with only the three follicular-occlusion phenotypes PSD, HS, ACNE_FG (just-identified)
#  M2a genetic multiple regression PSD ~ HS + BMI                               (is PSD-HS sharing independent of adiposity?)
#  M2b PSD ~ HS + ACNE_FG + BMI + WHRadjBMI + BALDP1                           (independent contributions of each trait)
# BALDP1 = UKB balding pattern 1 (no balding): a negative coefficient means a positive association with balding.
# Liability scale: PSD K=0.7%, HS K=1%, hospital acne K=1%; broad acne and BALDP1 use their sample prevalence.
# Outputs: results/genomicsem_rg_matrix.tsv, results/genomicsem_models.tsv, results/genomicsem_fit.tsv
suppressMessages({library(GenomicSEM); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
mu <- file.path(base, "data/processed/ldsc_munged")
ref <- paste0(file.path(base, "data/ldsc_ref/eur_w_ld_chr"), "/")
tr <- data.table(
  name  = c("PSD", "HS", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI"),
  file  = c("PSD_META_EUR", "HS_GCST90726609", "ACNE_FinnGen_R13", "ACNE_GCST90245818", "MPB_GCST90043616", "BMI_Pulit2019", "WHRadjBMI_Pulit2019"),
  sprev = c(6771 / 1381000, 3941 / 1439544, 4812 / 479708, 34422 / 399413, 66172 / 207036, NA, NA),
  pprev = c(0.007, 0.01, 0.01, 34422 / 399413, 66172 / 207036, NA, NA))
setwd(file.path(base, "results"))
LD <- ldsc(traits = file.path(mu, paste0(tr$file, ".sumstats.gz")), sample.prev = tr$sprev, population.prev = tr$pprev,
           ld = ref, wld = ref, trait.names = tr$name, stand = TRUE, ldsc.log = "genomicsem_ldsc")
rg <- as.data.table(LD$S_Stand, keep.rownames = FALSE); rg[, trait := tr$name]; setcolorder(rg, "trait")
fwrite(rg, file.path(base, "results/genomicsem_rg_matrix.tsv"), sep = "\t")
cat("\nGenomic SEM LDSC genetic correlation matrix:\n"); print(round(LD$S_Stand, 3))

pairs_cov <- function(v) paste(combn(v, 2, function(x) paste(x[1], "~~", x[2])), collapse = "\n")
models <- list(
  M1_common_factor_4 = "F =~ NA*PSD + HS + ACNE_FG + ACNE_TL\nF ~~ 1*F",
  M1b_common_factor_3 = "F =~ NA*PSD + HS + ACNE_FG\nF ~~ 1*F",
  M2a_PSD_on_HS_BMI = paste("PSD ~ HS + BMI", pairs_cov(c("HS", "BMI")), sep = "\n"),
  M2b_PSD_on_all = paste("PSD ~ HS + ACNE_FG + BMI + WHRadjBMI + BALDP1", pairs_cov(c("HS", "ACNE_FG", "BMI", "WHRadjBMI", "BALDP1")), sep = "\n"))
res <- list(); fit <- list()
for (m in names(models)) {
  um <- tryCatch(usermodel(LD, estimation = "DWLS", model = models[[m]], CFIcalc = TRUE, std.lv = FALSE, imp_cov = FALSE),
                 error = function(e) { cat("model", m, "failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(um)) next
  r <- as.data.table(um$results); r[, model := m]; res[[m]] <- r
  f <- as.data.table(um$modelfit); f[, model := m]; fit[[m]] <- f
  cat("\n====", m, "====\n"); print(um$results); print(um$modelfit)
}
fwrite(rbindlist(res, fill = TRUE), file.path(base, "results/genomicsem_models.tsv"), sep = "\t")
fwrite(rbindlist(fit, fill = TRUE), file.path(base, "results/genomicsem_fit.tsv"), sep = "\t")
cat("GENOMICSEM_DONE\n")
