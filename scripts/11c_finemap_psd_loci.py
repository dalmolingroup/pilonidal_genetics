#!/usr/bin/env python3
"""
Step 11c: Statistical fine-mapping of PSD meta-analysis loci.
  (1) Approximate Bayes Factor (Wakefield 2009) 95% credible sets on META_ALL and META_EUR (LD-free).
  (2) SuSiE-RSS (Zou et al. 2022; susieR) on META_EUR using 1000 Genomes Phase 3 EUR (n=503, GRCh38) LD from plink2,
      with the kriging_rss LD/z-score consistency diagnostic.
  (3) Functional annotation of credible-set variants: Ensembl VEP (consequence, regulatory feature, TF motif disruption)
      and GTEx v8 skin / subcutaneous adipose single-tissue eQTL lookup.
Inputs : results/psd_meta_{ALL,EUR}_loci.tsv, data/processed/meta/psd_meta_{ALL,EUR}.tsv.gz, data/raw/1kg_hg38/
Outputs: results/psd_finemap_credible_sets.tsv, results/psd_finemap_locus_summary.tsv, results/psd_finemap_vep_annotation.tsv,
         results/psd_finemap_gtex_eqtl.tsv, figures/psd_finemap_locuszoom_*.png
"""
import json
import math
import os
import subprocess
import time
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
from scipy.stats import norm

BASE = Path(__file__).resolve().parent.parent
META = BASE / "data" / "processed" / "meta"
RES = BASE / "results"
FIG = BASE / "figures"
KG = BASE / "data" / "raw" / "1kg_hg38"
FM = BASE / "data" / "processed" / "finemap"
FM.mkdir(parents=True, exist_ok=True)
PLINK2 = os.environ.get("PLINK2", "/tmp/claude-1000/-home-jvfe-dev-pesquisa-minor-contribs-pilonidal-genetics/5cd42c41-c3f4-4a05-a060-fe0ffae56124/scratchpad/bin/plink2")
WINDOW = 500_000
W_PRIOR = 0.2 ** 2          # prior variance of log-OR effect (sd = 0.2)
P_LOCUS = float(os.environ.get("FINEMAP_P", "5e-8"))   # loci to fine-map (default genome-wide significant)
N_EFF_EUR = 4 / (1 / 1433 + 1 / 448820) + 4 / (1 / 1466 + 1 / 456974) + 4 / (1 / 3872 + 1 / 468418)
GTEX_TISSUES = ["Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum"]


def http_json(url, payload=None, retries=3):
    for i in range(retries):
        try:
            data = json.dumps(payload).encode() if payload is not None else None
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == retries - 1:
                print(f"    HTTP failed {url[:80]}: {e}", flush=True)
                return None
            time.sleep(3 * (i + 1))


def abf_credible_set(df):
    """Wakefield ABF with prior W on log-OR; returns df with abf, pip, in_cs95."""
    v = df["se"].values ** 2
    z = df["beta"].values / df["se"].values
    log_abf = 0.5 * np.log(v / (v + W_PRIOR)) + 0.5 * z ** 2 * W_PRIOR / (v + W_PRIOR)   # log BF for H1 vs H0
    log_abf -= log_abf.max()
    pip = np.exp(log_abf) / np.exp(log_abf).sum()
    df = df.assign(abf_pip=pip).sort_values("abf_pip", ascending=False)
    df["abf_cum"] = df["abf_pip"].cumsum()
    cs = np.zeros(len(df), dtype=bool)
    idx = int(np.searchsorted(df["abf_cum"].values, 0.95)) + 1
    cs[:idx] = True
    df["abf_cs95"] = cs
    return df


def main():
    con = duckdb.connect()
    loci = {m: pd.read_csv(RES / f"psd_meta_{m}_loci.tsv", sep="\t") for m in ("ALL", "EUR")}
    # Union of GWS loci (positions) across both metas, de-duplicated by 500 kb
    cand = pd.concat([loci[m][loci[m].p_meta < P_LOCUS].assign(src=m) for m in loci]).sort_values("p_meta")
    regions = []
    for r in cand.itertuples():
        if any(x["chr"] == str(r.chr) and abs(x["pos"] - int(r.pos)) < WINDOW for x in regions):
            continue
        regions.append({"chr": str(r.chr), "pos": int(r.pos), "lead_rsid": r.rsid, "lead_gene": r.nearest_coding_gene, "lead_p_src": r.src, "lead_p": r.p_meta})
    print(f"[1/5] {len(regions)} loci to fine-map (p < {P_LOCUS})", flush=True)
    if not regions:
        print("No loci reach the fine-mapping threshold; nothing to do."); return

    # -------- prepare 1000G EUR plink2 subset restricted to all loci (one plink2 call) --------
    print("[2/5] Extracting 1000G EUR genotypes for all loci windows...", flush=True)
    pgen = KG / "all_hg38.pgen"
    if not pgen.exists():
        subprocess.run([PLINK2, "--zst-decompress", str(KG / "all_hg38.pgen.zst"), str(pgen)], check=True)
    psam = pd.read_csv(KG / "hg38_corrected.psam", sep="\t")
    eur = psam[psam["SuperPop"] == "EUR"]["#IID"]
    (FM / "eur.keep").write_text("\n".join(f"{i}" for i in eur) + "\n")
    with open(FM / "loci.ranges", "w") as f:
        for r in regions:
            f.write(f"{r['chr']}\t{max(1, r['pos'] - WINDOW)}\t{r['pos'] + WINDOW}\tlocus_{r['chr']}_{r['pos']}\n")
    sub = FM / "kg_eur_loci"
    if not (sub.with_suffix(".pgen")).exists():
        subprocess.run([PLINK2, "--pfile", str(KG / "all_hg38"), "vzs", "--psam", str(KG / "hg38_corrected.psam"), "--keep", str(FM / "eur.keep"),
                        "--extract", "range", str(FM / "loci.ranges"), "--snps-only", "just-acgt", "--max-alleles", "2", "--maf", "0.005", "--rm-dup", "force-first",
                        "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "30", "missing", "--make-pgen", "--out", str(sub), "--threads", "8"], check=True,
                       stdout=subprocess.DEVNULL)
    pvar = pd.read_csv(sub.with_suffix(".pvar"), sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 4], names=["chr", "pos", "id", "ref", "alt"], dtype=str)
    pvar["pos"] = pvar["pos"].astype(int)
    pvar["key"] = pvar.apply(lambda r: f"{r.chr}:{r.pos}:{min(r.ref, r.alt)}:{max(r.ref, r.alt)}", axis=1)
    print(f"   {len(pvar):,} 1000G EUR variants across loci windows", flush=True)

    meta = {m: con.execute(f"SELECT chr, CAST(pos AS BIGINT) pos, a1, a2, rsid, CAST(k AS INT) k, CAST(beta_meta AS DOUBLE) beta, CAST(se_meta AS DOUBLE) se, CAST(p_meta AS DOUBLE) p, CAST(eaf_meta AS DOUBLE) eaf, direction FROM read_csv('{META}/psd_meta_{m}.tsv.gz', delim='\\t', header=true, all_varchar=true, quote='', escape='', comment='', strict_mode=false, ignore_errors=true, max_line_size=1000000) WHERE CAST(k AS INT) >= 2").df() for m in ("ALL", "EUR")}

    cs_rows, summary_rows = [], []
    rscript = """
suppressMessages(library(susieR))
args <- commandArgs(trailingOnly=TRUE)
z <- scan(args[1], quiet=TRUE); R <- as.matrix(read.table(args[2])); n <- as.numeric(args[3]); out <- args[4]
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
    (FM / "run_susie.R").write_text(rscript)

    print("[3/5] ABF + SuSiE fine-mapping per locus...", flush=True)
    for reg in regions:
        c, p0 = reg["chr"], reg["pos"]
        lo, hi = p0 - WINDOW, p0 + WINDOW
        tag = f"chr{c}_{p0}_{reg['lead_gene'] or reg['lead_rsid']}"
        for m in ("ALL", "EUR"):
            d = meta[m][(meta[m].chr == c) & (meta[m].pos.between(lo, hi))].copy()
            if d.empty:
                continue
            d = abf_credible_set(d)
            d["meta"] = m; d["locus"] = tag
            d["susie_pip"] = np.nan; d["susie_cs"] = np.nan; d["susie_cs_purity"] = np.nan; d["susie_ld_mismatch_flag"] = np.nan
            if m == "EUR":
                # ---- SuSiE-RSS with 1000G EUR LD ----
                d["key"] = d.chr + ":" + d.pos.astype(str) + ":" + d.a1 + ":" + d.a2
                mm = d.merge(pvar[["key", "id", "ref", "alt"]], on="key", how="inner")
                mm = mm[mm.p.notna()].drop_duplicates("key")
                if len(mm) >= 20:
                    ids = FM / f"{tag}.ids"; ids.write_text("\n".join(mm.id) + "\n")
                    pref = FM / f"{tag}_ld"
                    subprocess.run([PLINK2, "--pfile", str(sub), "--extract", str(ids), "--r-unphased", "square", "ref-based", "--out", str(pref), "--threads", "8"],
                                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    ld_ids = [l.strip() for l in open(str(pref) + ".unphased.vcor1.vars")]
                    R = np.loadtxt(str(pref) + ".unphased.vcor1")
                    mm = mm.set_index("id").loc[ld_ids].reset_index()
                    # correlations are REF-based; align z to the 1000G REF allele then flip back for reporting
                    z = (mm.beta / mm.se).values * np.where(mm.a1 == mm.ref, 1.0, -1.0)
                    np.savetxt(FM / f"{tag}.z", z)
                    np.savetxt(FM / f"{tag}.R", R, fmt="%.5f")
                    out = FM / f"{tag}.susie.tsv"
                    try:
                        subprocess.run(["Rscript", str(FM / "run_susie.R"), str(FM / f"{tag}.z"), str(FM / f"{tag}.R"), str(int(N_EFF_EUR)), str(out)], check=True,
                                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=1800)
                        s = pd.read_csv(out, sep="\t")
                        s["key"] = mm.key.values
                        d = d.merge(s.rename(columns={"pip": "susie_pip_", "cs": "susie_cs_", "purity": "susie_cs_purity_", "ldflag": "susie_ld_mismatch_flag_", "s_rss": "susie_s_rss"}), on="key", how="left")
                        for col in ("susie_pip", "susie_cs", "susie_cs_purity", "susie_ld_mismatch_flag"):
                            d[col] = d[col + "_"]; d.drop(columns=[col + "_"], inplace=True)
                    except subprocess.CalledProcessError as e:
                        print(f"    SuSiE failed for {tag}: {e.stderr.decode()[-300:]}", flush=True)
                    except subprocess.TimeoutExpired:
                        print(f"    SuSiE timed out for {tag}", flush=True)
            n_cs = int(d.abf_cs95.sum())
            top = d.iloc[0]
            srow = {"locus": tag, "meta": m, "chr": c, "window_start": lo, "window_end": hi, "lead_rsid": reg["lead_rsid"], "lead_gene": reg["lead_gene"],
                    "n_variants_window": len(d), "abf_top_rsid": top.rsid, "abf_top_pos": int(top.pos), "abf_top_pip": round(float(top.abf_pip), 3), "abf_cs95_size": n_cs,
                    "abf_cs95_span_bp": int(d[d.abf_cs95].pos.max() - d[d.abf_cs95].pos.min()) if n_cs else 0}
            if m == "EUR" and d.susie_pip.notna().any():
                ncs = int(d.susie_cs.dropna().nunique())
                srow.update({"susie_n_cs": ncs, "susie_s_rss": float(d.get("susie_s_rss", pd.Series([np.nan])).dropna().iloc[0]) if "susie_s_rss" in d and d["susie_s_rss"].notna().any() else np.nan,
                             "susie_top_rsid": d.sort_values("susie_pip", ascending=False).iloc[0].rsid, "susie_top_pip": round(float(d.susie_pip.max()), 3),
                             "susie_cs_sizes": ";".join(str(int((d.susie_cs == i).sum())) for i in sorted(d.susie_cs.dropna().unique())),
                             "susie_n_ld_mismatch_flags": int(d.susie_ld_mismatch_flag.fillna(False).astype(bool).sum())})
            summary_rows.append(srow)
            keep = d[(d.abf_cs95) | (d.susie_cs.notna()) | (d.susie_pip > 0.05)]
            cs_rows.append(keep)
            print(f"   {tag} META_{m}: window n={len(d)}, ABF top {top.rsid} PIP={top.abf_pip:.2f}, CS95 size={n_cs}" + (f", SuSiE CS={srow.get('susie_n_cs')} s_rss={srow.get('susie_s_rss')}" if m == "EUR" else ""), flush=True)

    cs = pd.concat(cs_rows, ignore_index=True)
    cs_cols = ["locus", "meta", "chr", "pos", "rsid", "a1", "a2", "eaf", "k", "direction", "beta", "se", "p", "abf_pip", "abf_cs95", "susie_pip", "susie_cs", "susie_cs_purity", "susie_ld_mismatch_flag"]
    cs = cs[[c for c in cs_cols if c in cs.columns]].sort_values(["locus", "meta", "abf_pip"], ascending=[True, True, False])
    cs.to_csv(RES / "psd_finemap_credible_sets.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(summary_rows).to_csv(RES / "psd_finemap_locus_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)

    # -------- VEP annotation of credible-set variants (ALL-meta ABF CS95 or SuSiE CS or PIP>0.05) --------
    print("[4/5] Ensembl VEP annotation (regulatory + motif) of credible-set variants...", flush=True)
    ann_vars = cs[(cs.abf_cs95) | (cs.susie_cs.notna()) | (cs.abf_pip > 0.05)].drop_duplicates(["chr", "pos", "a1", "a2"])
    ann_vars = ann_vars.sort_values("abf_pip", ascending=False).groupby("locus").head(60)
    vep_rows = []
    def vep_region_string(r):
        # VEP region format: chr start end ref/alt strand ; we use a1 as ref proxy and a2 as alt (alleles are canonical-sorted, VEP tolerates allele mismatch by reporting both)
        return f"{r.chr} {r.pos} {r.pos + max(len(r.a1), len(r.a2)) - 1} {r.a1}/{r.a2} + {r.rsid or 'NA'}"
    recs = list(ann_vars.itertuples())
    for i in range(0, len(recs), 150):
        batch = recs[i:i + 150]
        res = http_json("https://rest.ensembl.org/vep/homo_sapiens/region?regulatory=1&canonical=1&symbol=1", {"variants": [vep_region_string(r) for r in batch]})
        if not res:
            continue
        for v in res:
            inp = v.get("input", ""); parts = inp.split()
            key = (parts[0], int(parts[1]), parts[3].split("/")[0], parts[3].split("/")[1]) if len(parts) >= 4 else None
            tcs = v.get("transcript_consequences", []) or []
            canon = [t for t in tcs if t.get("canonical")] or tcs
            genes = ",".join(dict.fromkeys(t.get("gene_symbol", "") for t in canon if t.get("gene_symbol")))
            regs = v.get("regulatory_feature_consequences", []) or []
            motifs = v.get("motif_feature_consequences", []) or []
            vep_rows.append({"chr": key[0] if key else None, "pos": key[1] if key else None, "a1": key[2] if key else None, "a2": key[3] if key else None,
                             "most_severe_consequence": v.get("most_severe_consequence"), "genes_canonical": genes,
                             "distance_to_gene": ";".join(str(t.get("distance", 0)) for t in canon[:3]),
                             "regulatory_features": ";".join(f"{r.get('regulatory_feature_id')}:{r.get('biotype')}" for r in regs),
                             "n_motif_features": len(motifs),
                             "tf_motifs": ";".join(f"{m.get('transcription_factors', [''])[0] if m.get('transcription_factors') else m.get('motif_name','')}(score_change={m.get('motif_score_change')}{',HIGH_INFO' if m.get('high_inf_pos')=='Y' else ''})" for m in motifs[:8]),
                             "max_abs_motif_score_change": max([abs(float(m.get("motif_score_change", 0) or 0)) for m in motifs], default=0.0),
                             "high_inf_pos_motif": any(m.get("high_inf_pos") == "Y" for m in motifs)})
        time.sleep(1)
    vep = pd.DataFrame(vep_rows)
    if len(vep):
        vep = ann_vars[["locus", "meta", "chr", "pos", "rsid", "a1", "a2", "p", "abf_pip", "susie_pip", "susie_cs"]].merge(vep, on=["chr", "pos", "a1", "a2"], how="left")
        vep.to_csv(RES / "psd_finemap_vep_annotation.tsv", sep="\t", index=False, float_format="%.4g")
        print(f"   annotated {len(vep)} credible-set variants; with TF-motif disruption: {int((vep.n_motif_features.fillna(0) > 0).sum())}; in regulatory features: {int((vep.regulatory_features.fillna('') != '').sum())}", flush=True)

    # -------- GTEx v8 single-tissue eQTL lookup for top credible-set variants --------
    print("[5/5] GTEx v8 eQTL lookup (skin / subcutaneous adipose) for top credible-set variants...", flush=True)
    gtex_rows = []
    top_vars = ann_vars.sort_values("abf_pip", ascending=False).groupby("locus").head(15)
    for r in top_vars.itertuples():
        for ref, alt in ((r.a1, r.a2), (r.a2, r.a1)):
            vid = f"chr{r.chr}_{r.pos}_{ref}_{alt}_b38"
            res = http_json(f"https://gtexportal.org/api/v2/association/singleTissueEqtl?variantId={vid}&datasetId=gtex_v8&itemsPerPage=250")
            data = (res or {}).get("data", []) if isinstance(res, dict) else []
            if data:
                for e in data:
                    if e.get("tissueSiteDetailId") in GTEX_TISSUES:
                        gtex_rows.append({"locus": r.locus, "rsid": r.rsid, "chr": r.chr, "pos": r.pos, "variantId": vid, "abf_pip": r.abf_pip, "susie_pip": r.susie_pip,
                                          "gene": e.get("geneSymbol"), "gencodeId": e.get("gencodeId"), "tissue": e.get("tissueSiteDetailId"), "nes": e.get("nes"), "pValue": e.get("pValue")})
                break
            time.sleep(0.3)
    gtex = pd.DataFrame(gtex_rows)
    gtex.to_csv(RES / "psd_finemap_gtex_eqtl.tsv", sep="\t", index=False, float_format="%.4g")
    if len(gtex):
        print(gtex.sort_values("pValue").groupby(["locus", "gene", "tissue"]).head(1).head(40).to_string(index=False), flush=True)
    print("FINEMAP_DONE", flush=True)


if __name__ == "__main__":
    main()
