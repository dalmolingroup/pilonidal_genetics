#!/usr/bin/env python3
"""
Step 01f: SuSiE-RSS sensitivity analysis of the META_EUR loci with UK Biobank LD (Weissbrod et al. 2020; 337k unrelated
British-ancestry participants, GRCh37, 3 Mb windows from s3://broad-alkesgroup-ukbb-ld/UKBB_LD) instead of the
1000 Genomes EUR LD (n=503) used in step 01c. SuSiE settings, windows (GRCh38 lead +/- 500 kb) and N are identical to 01c.
  - meta variants are matched to the UKB LD variants by rsID and allele pair; z-scores are aligned to the UKB allele1
  - GRCh37 window centres come from the median GRCh37-GRCh38 offset of window variants in the LAVA UKB .info files
  - allele-alignment check: share of variants in strong LD with the lead (|r| > 0.8) whose z sign agrees with r * z_lead
Inputs : results/psd_finemap_locus_summary.tsv (01c), data/processed/meta/psd_meta_EUR.tsv.gz, data/lava/ref/*.info
Outputs: results/psd_finemap_susie_ukbld.tsv (variant level), results/psd_finemap_susie_ukbld_summary.tsv
         (per locus, with the 1000G SuSiE columns from 01c alongside); LD files cached in data/raw/ukbb_ld/
"""
import subprocess
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
from scipy import sparse

BASE = Path(__file__).resolve().parent.parent
META = BASE / "data" / "processed" / "meta"
RES = BASE / "results"
LAVA_REF = BASE / "data" / "lava" / "ref"
LD = BASE / "data" / "raw" / "ukbb_ld"
FM = BASE / "data" / "processed" / "finemap_ukbld"
LD.mkdir(parents=True, exist_ok=True)
FM.mkdir(parents=True, exist_ok=True)
URL = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD"
N_EFF_EUR = 4 / (1 / 1433 + 1 / 448820) + 4 / (1 / 1466 + 1 / 456974) + 4 / (1 / 3872 + 1 / 468418)   # as in 01c

RSCRIPT = """
suppressMessages(library(susieR))
args <- commandArgs(trailingOnly=TRUE)
z <- scan(args[1], quiet=TRUE); p <- length(z)
con <- file(args[2], "rb"); R <- matrix(readBin(con, "double", p * p, size=8, endian="little"), p, p); close(con)
n <- as.numeric(args[3]); out <- args[4]
R[is.na(R)] <- 0; diag(R) <- 1
fit <- tryCatch(susie_rss(z=z, R=R, n=n, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=TRUE, max_iter=500, verbose=FALSE), error=function(e) NULL)
if (is.null(fit)) fit <- susie_rss(z=z, R=R, n=n, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=FALSE, max_iter=500, verbose=FALSE)
pip <- susie_get_pip(fit); cs <- rep(NA_integer_, length(z))
if (!is.null(fit$sets$cs)) for (i in seq_along(fit$sets$cs)) cs[fit$sets$cs[[i]]] <- i
purity <- rep(NA_real_, length(z)); if (!is.null(fit$sets$purity)) for (i in seq_along(fit$sets$cs)) purity[fit$sets$cs[[i]]] <- fit$sets$purity$min.abs.corr[i]
cond <- tryCatch(kriging_rss(z, R, n=n), error=function(e) NULL)
lam <- tryCatch(estimate_s_rss(z, R, n=n), error=function(e) NA)
flag <- if (!is.null(cond)) cond$conditional_dist$logLR > 2 & abs(z) > 2 else rep(NA, length(z))
write.table(data.frame(pip=pip, cs=cs, purity=purity, ldflag=flag, s_rss=lam), out, sep="\\t", quote=FALSE, row.names=FALSE)
"""


def url_exists(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            return r.status == 200
    except Exception:
        return False


def fetch(name):
    dest = LD / name
    if not (dest.exists() and (LD / f"{name}.done").exists()):
        subprocess.run(["curl", "-sSL", "--retry", "6", "--retry-delay", "15", "-C", "-", "-o", str(dest), f"{URL}/{name}"], check=True)
        (LD / f"{name}.done").touch()
    return dest


def ld_window(chrom, centre37):
    """3 Mb UKB LD file (starts at 1 + k Mb) whose centre is closest to centre37; falls back to neighbours at chromosome ends."""
    k0 = max(0, int(round((centre37 - 1_500_000) / 1_000_000)))
    for k in (k0, k0 - 1, k0 + 1, k0 - 2, k0 + 2):
        if k < 0:
            continue
        prefix = f"chr{chrom}_{k * 1_000_000 + 1}_{k * 1_000_000 + 3_000_001}"
        if url_exists(f"{URL}/{prefix}.gz"):
            return prefix
    return None


def main():
    con = duckdb.connect()
    loci = pd.read_csv(RES / "psd_finemap_locus_summary.tsv", sep="\t", dtype={"chr": str})
    loci = loci[loci.meta == "EUR"]
    loci = loci[~((loci.chr == "6") & (loci.window_end > 25_000_000) & (loci.window_start < 34_000_000))]   # MHC excluded
    print(f"[1/3] {len(loci)} non-MHC loci", flush=True)
    meta = con.execute(f"SELECT chr, CAST(pos AS BIGINT) pos, a1, a2, rsid, CAST(beta_meta AS DOUBLE) beta, CAST(se_meta AS DOUBLE) se, CAST(p_meta AS DOUBLE) p FROM read_csv('{META}/psd_meta_EUR.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) WHERE CAST(k AS INT) >= 2").df()
    (FM / "run_susie.R").write_text(RSCRIPT)

    var_rows, summary_rows = [], []
    lava_info = {}
    for reg in loci.itertuples():
        c, lo, hi, tag = reg.chr, int(reg.window_start), int(reg.window_end), reg.locus
        d = meta[(meta.chr == c) & meta.pos.between(lo, hi) & meta.p.notna() & meta.rsid.str.startswith("rs", na=False)].drop_duplicates("rsid")
        if c not in lava_info:
            lava_info[c] = pd.read_csv(LAVA_REF / f"lava-ukb-v1.1_chr{c}.info", sep="\t", usecols=["SNP", "POS"])
        off = d.merge(lava_info[c], left_on="rsid", right_on="SNP")
        if off.empty:
            print(f"   {tag}: no GRCh37 anchor variants, skipped", flush=True); continue
        centre37 = int((lo + hi) / 2 + np.median(off.POS - off.pos))
        prefix = ld_window(c, centre37)
        if prefix is None:
            print(f"   {tag}: no UKB LD file found near chr{c}:{centre37}", flush=True); continue
        snps = pd.read_csv(fetch(f"{prefix}.gz"), sep="\t", dtype={"chromosome": str})
        snps["idx"] = np.arange(len(snps))
        mm = d.merge(snps, on="rsid")
        mm = mm[((mm.a1 == mm.allele1) & (mm.a2 == mm.allele2)) | ((mm.a1 == mm.allele2) & (mm.a2 == mm.allele1))].drop_duplicates("rsid")
        if len(mm) < 20:
            print(f"   {tag}: only {len(mm)} matched variants, skipped", flush=True); continue
        mm = mm.sort_values("idx").reset_index(drop=True)
        U = sparse.load_npz(fetch(f"{prefix}.npz")).tocsr()
        assert sparse.tril(U, k=-1).nnz == 0 or sparse.triu(U, k=1).nnz == 0, "expected a triangular LD matrix"
        sub = U[mm.idx.values][:, mm.idx.values].toarray()
        R = sub + sub.T
        np.fill_diagonal(R, 1.0)
        # z for the meta a1 allele, aligned to the UKB allele1 so that z and R refer to the same allele at every variant
        z = (mm.beta / mm.se).values * np.where(mm.a1 == mm.allele1, 1.0, -1.0)
        lead = int(np.argmax(np.abs(z)))
        strong = (np.abs(R[lead]) > 0.8) & (np.arange(len(z)) != lead)
        sign_ok = float(np.mean(np.sign(z[strong]) == np.sign(R[lead, strong] * z[lead]))) if strong.any() else np.nan
        np.savetxt(FM / f"{tag}.z", z)
        R.astype("<f8").tofile(FM / f"{tag}.R.bin")
        out = FM / f"{tag}.susie.tsv"
        try:
            subprocess.run(["Rscript", str(FM / "run_susie.R"), str(FM / f"{tag}.z"), str(FM / f"{tag}.R.bin"), str(int(N_EFF_EUR)), str(out)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=3600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f"   {tag}: SuSiE failed ({getattr(e, 'stderr', b'')[-300:]})", flush=True); continue
        (FM / f"{tag}.R.bin").unlink()
        s = pd.concat([mm[["chr", "pos", "rsid", "a1", "a2", "beta", "se", "p", "position"]].rename(columns={"position": "pos_b37"}), pd.read_csv(out, sep="\t")], axis=1)
        s.insert(0, "locus", tag)
        var_rows.append(s[(s.cs.notna()) | (s.pip > 0.01)])
        top = s.sort_values("pip", ascending=False).iloc[0]
        summary_rows.append({"locus": tag, "lead_rsid": reg.lead_rsid, "ld_file": prefix, "n_window": len(d), "n_matched": len(mm), "lead_sign_concordance": round(sign_ok, 3),
                             "ukb_susie_n_cs": int(s.cs.dropna().nunique()), "ukb_susie_cs_sizes": ";".join(str(int((s.cs == i).sum())) for i in sorted(s.cs.dropna().unique())),
                             "ukb_susie_top_rsid": top.rsid, "ukb_susie_top_pip": round(float(top.pip), 3), "ukb_s_rss": float(s.s_rss.iloc[0]),
                             "ukb_n_ld_mismatch_flags": int(s.ldflag.fillna(False).astype(bool).sum()),
                             "kg_susie_n_cs": reg.susie_n_cs, "kg_susie_cs_sizes": reg.susie_cs_sizes, "kg_susie_top_rsid": reg.susie_top_rsid, "kg_susie_top_pip": reg.susie_top_pip,
                             "kg_s_rss": reg.susie_s_rss, "kg_n_ld_mismatch_flags": reg.susie_n_ld_mismatch_flags})
        r = summary_rows[-1]
        print(f"   {tag}: matched {len(mm)}/{len(d)}, sign concordance {sign_ok:.2f}, CS {r['ukb_susie_n_cs']} (1000G {reg.susie_n_cs}), top {top.rsid} PIP {top.pip:.2f} (1000G {reg.susie_top_rsid} {reg.susie_top_pip}), s {r['ukb_s_rss']:.3f} (1000G {reg.susie_s_rss})", flush=True)

    print("[2/3] writing results", flush=True)
    pd.concat(var_rows, ignore_index=True).to_csv(RES / "psd_finemap_susie_ukbld.tsv", sep="\t", index=False, float_format="%.4g")
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(RES / "psd_finemap_susie_ukbld_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print("[3/3] summary", flush=True)
    print(summ.drop(columns=["ld_file"]).to_string(index=False), flush=True)
    print("FINEMAP_UKBLD_DONE", flush=True)


if __name__ == "__main__":
    main()
