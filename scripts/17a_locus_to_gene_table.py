#!/usr/bin/env python3
"""Step 17a: locus-to-gene evidence table for the 25 non-MHC PSD meta-analysis loci.

Collects, per locus, the evidence already produced by steps 11-15 and scores candidate genes transparently:
  nearest gene (+1), top MAGMA gene among genes within 250 kb if p < 1e-3 (+1), coding variant in the fine-mapped credible set
  (PIP >= 0.01; +2 if PIP >= 0.1, else +1), GTEx eQTL colocalization (+2 if PP.H4 >= 0.8, +1 if 0.5-0.8), Open Targets L2G top gene (+1).
A protein-coding gene with functional evidence (coding PIP >= 0.1 or eQTL PP.H4 >= 0.8) outranks proximity.
Locus-level context (not gene-specific): coloc with HS / acne / balding / adiposity (PP.H4 >= 0.8, META_EUR) and
LAVA local rg (FDR < 0.05; GRCh37 blocks, lead positions looked up via the Ensembl GRCh37 REST API and cached).
Inputs : results/psd_meta_locus_classification.tsv, psd_meta_loci_intercept_corrected.tsv, psd_finemap_locus_summary.tsv,
         psd_finemap_vep_by_rsid.tsv, magma_genes_PSD_META_EUR.tsv, eqtl_coloc_psd_gtex.tsv, coloc_psd_vs_traits.tsv,
         lava_bivariate_psd_noMHC.tsv, opentargets_psd_known_credible_sets.tsv
Output : results/psd_locus_to_gene.tsv
"""
import json, urllib.request
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
R = BASE / "results"
rd = lambda f: pd.read_csv(R / f, sep="\t")
CODING = {"missense_variant", "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost", "stop_lost"}
# Ensembl GRCh38 (checked 2026-10-21): only TANAR (HGNC:56695) and RMDN2-AS1 (HGNC:41150) have symbols; the others are unnamed novel transcripts
LNC = {"ENSG00000287727": "ENSG00000287727, antisense to PEX14 (lncRNA)", "ENSG00000223991": "TANAR (lncRNA)",
       "ENSG00000236838": "ENSG00000236838, antisense to SMG6 (lncRNA)", "ENSG00000232821": "ENSG00000232821, overlaps TWIST1 (lncRNA)",
       "ENSG00000235848": "RMDN2-AS1 (lncRNA)"}
TISSUE = {"Skin_sun_exposed": "skin-SE", "Skin_not_sun_exposed": "skin-NSE", "Adipose_subcutaneous": "adipose-SC",
          "Adipose_visceral": "adipose-VO", "Fibroblast_cultured": "fibroblast"}

loci = rd("psd_meta_locus_classification.tsv")
loci = loci[loci.is_mhc.astype(str).str.lower() == "false"].reset_index(drop=True)
icorr = rd("psd_meta_loci_intercept_corrected.tsv")[["chr", "pos", "p_intercept_corrected", "gws_after_correction"]]
loci = loci.merge(icorr, on=["chr", "pos"], how="left").rename(columns={"or": "odds_ratio"})
fm = rd("psd_finemap_locus_summary.tsv"); fm["pos"] = fm.locus.str.split("_").str[1].astype(int)
vep = rd("psd_finemap_vep_by_rsid.tsv"); vep["pos_l"] = vep.locus.str.split("_").str[1].astype(int)
magma = rd("magma_genes_PSD_META_EUR.tsv").sort_values("P").drop_duplicates("SYMBOL").set_index("SYMBOL")
eq = rd("eqtl_coloc_psd_gtex.tsv")
eq["gene"] = [LNC.get(i, n if isinstance(n, str) and n else i) for i, n in zip(eq.gene_id, eq.gene_name)]
co = rd("coloc_psd_vs_traits.tsv"); co = co[(co.psd_set == "PSD_META_EUR") & (co.region_source == "PSD")]
lava = rd("lava_bivariate_psd_noMHC.tsv"); lava = lava[lava.fdr < 0.05]
ot = rd("opentargets_psd_known_credible_sets.tsv")
# the chr7:18874234 indel lead has no rsID in the meta-analysis; take it from the Open Targets credible set at the same position
ot_rs = ot.set_index(["chr", "pos_b38"]).rsid
loci["rsid"] = [r if isinstance(r, str) else ot_rs.get((c, p)) for r, c, p in zip(loci.rsid, loci.chr, loci.pos)]

# GRCh37 positions of the lead variants (LAVA blocks are GRCh37)
cache = BASE / "data/processed/psd_leads_grch37.json"
if cache.exists():
    b37 = json.loads(cache.read_text())
else:
    req = urllib.request.Request("https://grch37.rest.ensembl.org/variation/homo_sapiens", method="POST",
                                 data=json.dumps({"ids": [r for r in loci.rsid if isinstance(r, str)]}).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    js = json.load(urllib.request.urlopen(req, timeout=60))
    b37 = {k: next((m["start"] for m in v["mappings"] if m["assembly_name"] == "GRCh37" and "_" not in m["seq_region_name"]), None)
           for k, v in js.items()}
    cache.write_text(json.dumps(b37))

# nearest protein-coding gene (GRCh38, +/-1 Mb) for gene-desert loci without a nearest gene in step 11d
cache38 = BASE / "data/processed/psd_desert_nearest_genes.json"
desert = {f"{c}:{p}": None for c, p, g in zip(loci.chr, loci.pos, loci.nearest_gene) if not isinstance(g, str)}
near38 = json.loads(cache38.read_text()) if cache38.exists() else {}
for k in desert:
    if k in near38: continue
    c, p = map(int, k.split(":"))
    url = f"https://rest.ensembl.org/overlap/region/human/{c}:{max(p - 1_000_000, 1)}-{p + 1_000_000}?feature=gene;biotype=protein_coding"
    gs = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60))
    d = sorted((max(g["start"] - p, p - g["end"], 0), g.get("external_name") or g["id"]) for g in gs)
    near38[k] = [d[0][1], d[0][0]] if d else None
cache38.write_text(json.dumps(near38))

def is_lnc(g):
    return "lncRNA" in g or g.startswith(("ENSG", "LINC")) or "-AS" in g

rows = []
for L in loci.itertuples():
    genes250 = [g for g in str(L.genes_250kb).split(",") if g and g != "nan"]
    score, why, functional = {}, {}, set()
    def add(g, pts, reason):
        score[g] = score.get(g, 0) + pts; why.setdefault(g, []).append(reason)
    nearest, dist = L.nearest_gene, 0
    if not isinstance(nearest, str):
        nearest, dist = near38.get(f"{L.chr}:{L.pos}") or (None, None)
    if nearest: add(nearest, 1, "nearest" if not dist else f"nearest ({dist / 1e3:.0f} kb)")
    mg = magma.reindex(genes250).dropna(subset=["P"])
    magma_top = mg.P.idxmin() if len(mg) else None
    if magma_top and mg.P.min() < 1e-3: add(magma_top, 1, f"MAGMA p={mg.P.min():.1g}")
    # fine-mapping (ABF on META_ALL, SuSiE on META_EUR); coding variants with PIP >= 0.01 in the credible sets
    f = fm[(fm.chr == L.chr) & (fm.pos == L.pos)].set_index("meta")
    fa, fe = (f.loc[m] if m in f.index else None for m in ("ALL", "EUR"))
    v = vep[(vep.chr == L.chr) & (vep.pos_l == L.pos)].sort_values("pip_max", ascending=False)
    top = v.iloc[0] if len(v) else None
    cod = v[v.most_severe_consequence.isin(CODING) & (v.pip_max >= 0.01)].drop_duplicates("rsid")
    coding = []
    for c in cod.itertuples():
        g = str(c.genes_canonical).split(",")[0]
        aa = str(c.aa_change).split(";")[0].split(":")[-1] if isinstance(c.aa_change, str) else ""
        coding.append(f"{g} {c.most_severe_consequence.replace('_variant', '')} {aa} {c.rsid} (PIP {c.pip_max:.2f})".replace("  ", " "))
        add(g, 2 if c.pip_max >= 0.1 else 1, f"coding PIP {c.pip_max:.2f}")
        if c.pip_max >= 0.1: functional.add(g)
    # eQTL colocalization (best PP.H4 per gene across tissues and both meta-analyses)
    e = eq[(eq.chr == L.chr) & (eq.pos == L.pos) & (eq["PP.H4"] >= 0.5)]
    eqtl_txt = []
    for g, d in e.groupby("gene"):
        h4 = d["PP.H4"].max(); best = d.loc[d["PP.H4"].idxmax()]
        tis = ",".join(sorted({TISSUE[t] for t in d[d["PP.H4"] >= min(0.8, h4)].tissue}))
        eqtl_txt.append((h4, f"{g} {h4:.2f} [{tis}; risk allele {best.risk_allele_expression}]"))
        add(g, 2 if h4 >= 0.8 else 1, f"eQTL coloc {h4:.2f}")
        if h4 >= 0.8: functional.add(g)
    eqtl_txt = [t for _, t in sorted(eqtl_txt, reverse=True)]
    o = ot[(ot.chr == L.chr) & ((ot.pos_b38 - L.pos).abs() < 5e5)]
    l2g = list(dict.fromkeys(s.split(",")[0] for s in o.l2g_top3.dropna()))
    for g in l2g: add(g, 1, "OT L2G")
    # locus-level context
    c = co[(co.chr == L.chr) & (co.pos == L.pos) & (co["PP.H4"] >= 0.8)].sort_values("PP.H4", ascending=False)
    coloc_txt = "; ".join(f"{t} {h:.2f}" for t, h in zip(c.trait, c["PP.H4"]))
    p37 = b37.get(L.rsid)
    lv = lava[(lava.chr == L.chr) & (lava.start <= p37) & (lava.stop >= p37)] if p37 else lava.iloc[0:0]
    lava_txt = "; ".join(f"{t} {r:+.2f}" for t, r in zip(lv.phen1, lv.rho))
    # prioritisation: a protein-coding gene with functional evidence (coding PIP >= 0.1 or eQTL PP.H4 >= 0.8) outranks
    # proximity-based scores; otherwise the highest evidence score wins (ties reported jointly)
    fpc = [g for g in functional if not is_lnc(g)]
    pool = fpc or list(score)
    best = max((score[g] for g in pool), default=0)
    prio = sorted(g for g in pool if score[g] == best) or ["none"]
    if fpc: tier = "functional (protein-coding)"
    elif functional: tier = "functional (lncRNA only)"
    else: tier = "supportive (2+ lines)" if best >= 2 else "single line"
    rows.append(dict(
        locus=f"chr{L.chr}:{L.pos}", lead_rsid=L.rsid, effect_allele=L.a1, OR=L.odds_ratio, p_META_ALL=L.p_META_ALL,
        p_META_EUR=L.p_META_EUR, status=L.status, gws_after_intercept_correction=L.gws_after_correction, prior_art=L.prior_art,
        nearest_gene=nearest, nearest_gene_distance_kb=round(dist / 1e3) if dist else 0, genes_250kb=",".join(genes250),
        finemap_top_ALL=None if fa is None else f"{fa.abf_top_rsid} (PIP {fa.abf_top_pip:.2f}; CS95 {fa.abf_cs95_size:.0f})",
        finemap_top_EUR_susie=None if fe is None or pd.isna(fe.susie_top_rsid) else f"{fe.susie_top_rsid} (PIP {fe.susie_top_pip:.2f})",
        top_variant_consequence=None if top is None else f"{top.rsid} {top.most_severe_consequence} {top.genes_canonical if isinstance(top.genes_canonical, str) else ''}".strip(),
        coding_variants_in_cs="; ".join(coding), magma_top_gene=magma_top, eqtl_coloc="; ".join(eqtl_txt), ot_l2g=",".join(l2g),
        prioritised_gene="/".join(prio), functional_candidates=",".join(sorted(functional)), evidence_score=best, locus_evidence_tier=tier,
        evidence="; ".join(f"{g}: {', '.join(why[g])}" for g in sorted(score, key=score.get, reverse=True)[:4]),
        trait_coloc_pph4_08=coloc_txt, lava_local_rg_fdr05=lava_txt))
out = pd.DataFrame(rows).sort_values("p_META_ALL")
out.to_csv(R / "psd_locus_to_gene.tsv", sep="\t", index=False)
pd.set_option("display.width", 300, "display.max_colwidth", 70)
print(out[["locus", "status", "nearest_gene", "prioritised_gene", "locus_evidence_tier", "evidence", "trait_coloc_pph4_08", "lava_local_rg_fdr05"]].to_string(index=False))
print(out.locus_evidence_tier.value_counts().to_string())
