#!/usr/bin/env python3
"""Tables for the manuscript, from results/: Table 1 ($PSD_PAPER/tables/table1.tex, \\input by manuscript.tex and
the submission copy Table1.xlsx) and the supplementary tables workbook (supplementary_tables.xlsx). Run locally."""
import os
from pathlib import Path
import numpy as np, pandas as pd

BASE = Path(__file__).resolve().parents[1]
PAPER = Path(os.environ.get("PSD_PAPER", BASE / "paper"))  # manuscript directory (kept off the public repo)
RES, OUT = Path(os.environ.get("PSD_RESULTS", BASE / "results")), PAPER / "tables"
OUT.mkdir(parents=True, exist_ok=True)

def fmt_p(p):
    m, e = f"{p:.1e}".split("e"); return f"{m}×10\\textsuperscript{{{int(e)}}}".replace("{-", "{−")

def it(g):
    return f"\\textit{{{g}}}"

l2g = pd.read_csv(RES / "psd_locus_to_gene.tsv", sep="\t")
sup = pd.read_csv(RES / "psd_lead_variants_supplementary.tsv", sep="\t")
ic = pd.read_csv(RES / "psd_meta_loci_intercept_corrected.tsv", sep="\t")
suggestive = set(ic.loc[~ic.gws_after_correction.astype(bool), "rsid"])
eur, allm = sup[sup.meta == "META_EUR"].set_index("locus"), sup[sup.meta == "META_ALL"].set_index("locus")

def gene(p, n):
    p = str(p)
    if "antisense" in p or p.startswith("ENSG"):
        return f"{it(n)} (lncRNA)"
    return "/".join(it(g.replace(" (lncRNA)", "")) + (" (lncRNA)" if "lncRNA" in g else "") for g in p.split("/"))

rows, plain = [], []
loc = l2g.locus.str.extract(r"chr(\d+):(\d+)").astype(int)
for i in np.lexsort((loc[1].values, loc[0].values)):
    r = l2g.iloc[i]; e, a = eur.loc[r.locus], allm.loc[r.locus]
    mark = ("$^{\\dagger}$" if r.lead_rsid in suggestive else "") + ("$^{\\ddagger}$" if e.k_cohorts < 3 else "")
    ci = lambda s: s
    rows.append({"Locus (GRCh38)": f"{loc.iloc[i, 0]}:{loc.iloc[i, 1]:,}", "Lead variant": f"{r.lead_rsid}{mark}",
                 "EA/OA": f"{e.effect_allele}/{e.other_allele}", "EAF": f"{e.eaf:.2f}",
                 "OR (95\\% CI), EUR": ci(e.or_fixed_95ci), "\\textit{p}, EUR": fmt_p(e.p_fixed),
                 "OR (95\\% CI), ALL": ci(a.or_fixed_95ci), "\\textit{p}, ALL": fmt_p(a.p_fixed),
                 "Status": "known" if r.status == "KNOWN" else "novel", "Nearest gene": it(r.nearest_gene),
                 "Prioritised gene": gene(r.prioritised_gene, r.nearest_gene)})
    plain.append({"Locus (GRCh38)": f"{loc.iloc[i, 0]}:{loc.iloc[i, 1]}", "Lead variant": r.lead_rsid,
                  "Footnote": ("†" if r.lead_rsid in suggestive else "") + ("‡" if e.k_cohorts < 3 else ""),
                  "Effect allele": e.effect_allele, "Other allele": e.other_allele, "EAF (EUR)": round(e.eaf, 3),
                  "OR (95% CI), EUR": e.or_fixed_95ci, "P, EUR": float(f"{e.p_fixed:.2g}"),
                  "OR (95% CI), ALL": a.or_fixed_95ci, "P, ALL": float(f"{a.p_fixed:.2g}"),
                  "Status": "known" if r.status == "KNOWN" else "novel", "Nearest gene": r.nearest_gene,
                  "Prioritised gene": str(r.prioritised_gene)})
t1 = pd.DataFrame(rows)
head = " & ".join(t1.columns) + r" \\"
body = "\n".join(" & ".join(map(str, r)) + r" \\" for r in t1.itertuples(index=False))
note = (r"EA/OA, effect/other allele. EAF, effect-allele frequency in META\_EUR. OR, odds ratio per effect allele from "
        r"fixed-effect inverse-variance meta-analysis (METAL). EUR, European meta-analysis (MVP European, UK Biobank, FinnGen). ALL, "
        r"multi-ancestry meta-analysis (MVP multi-ancestry, UK Biobank, FinnGen, Genes \& Health). Prioritised gene: "
        r"locus-to-gene scoring (Materials and Methods). $^{\dagger}$\,$5\times10^{-8} < p < 10^{-7}$ after LDSC-intercept "
        r"correction. $^{\ddagger}$\,Lead variant measured in two of the three European cohorts.")
tex = "\n".join([
    r"\begin{table}[p]", r"\centering", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
    r"\caption{Twenty-five genome-wide significant pilonidal sinus disease loci outside the MHC.}\label{tab:loci}",
    r"\begin{tabular}{llllrrrrlll}", r"\toprule", head, r"\midrule", body, r"\bottomrule", r"\end{tabular}",
    r"\par\smallskip\parbox{\linewidth}{\footnotesize " + note + "}", r"\end{table}", ""])
(OUT / "table1.tex").write_text(tex)
(OUT / "table1.md").unlink(missing_ok=True)
# EJHG wants each table as its own Excel workbook, with its title and footnotes
with pd.ExcelWriter(OUT / "Table1.xlsx") as xw:
    title = "Table 1. Twenty-five genome-wide significant pilonidal sinus disease loci outside the MHC."
    pd.DataFrame([[title]]).to_excel(xw, sheet_name="Table 1", index=False, header=False)
    pd.DataFrame(plain).to_excel(xw, sheet_name="Table 1", index=False, startrow=2)
    notes = ("EAF, effect-allele frequency in META_EUR. OR, odds ratio per effect allele from fixed-effect inverse-variance "
             "meta-analysis (METAL). EUR, MVP European + UK Biobank + FinnGen. ALL, MVP multi-ancestry + UK Biobank + FinnGen + "
             "Genes & Health. † 5×10-8 < P < 1×10-7 after LDSC-intercept correction. ‡ Lead variant measured in two of the "
             "three European cohorts.")
    pd.DataFrame([[notes]]).to_excel(xw, sheet_name="Table 1", index=False, header=False, startrow=len(plain) + 4)

# variant nomenclature for the EJHG "Related Manuscript File" (curated HGVS in data/psd_variant_hgvs.tsv)
hg = pd.read_csv(BASE / "data/psd_variant_hgvs.tsv", sep="\t", dtype=str, keep_default_na=False)
hg = hg.assign(c=hg.chr.astype(int), p=hg.pos_grch38.astype(int)).sort_values(["c", "p"])
esc = lambda x: x.replace("_", "\\_").replace(">", "\\textgreater{}")
lines = [" & ".join([r.rsid, esc(r.hgvs_g), esc(", ".join(x for x in (r.hgvs_c, r.hgvs_p) if x)) or "--",
                     esc(r.consequence), it(r.gene) if r.gene else "--", r.role]) + r" \\" for r in hg.itertuples()]
(OUT / "nomenclature.tex").write_text("\n".join([
    r"\begin{tabular}{@{}llllll@{}}", r"\toprule",
    r"dbSNP & Genomic HGVS (GRCh38) & Transcript / protein HGVS & Consequence (VEP) & Gene (VEP) & Role \\",
    r"\midrule", *lines, r"\bottomrule", r"\end{tabular}", ""]))

# S4: both robustness checks of the lead variants (LDSC-intercept correction and the cohorts other than FinnGen) per locus
nonfg = pd.read_csv(RES / "psd_meta_loci_nonfinngen_replication.tsv", sep="\t")
robust = ic.assign(locus="chr" + ic.chr.astype(str) + ":" + ic.pos.astype(str)).merge(
    nonfg.drop(columns=["status", "gene", "rsid"]).rename(columns={"meta": "nonFG_meta", "p_meta": "nonFG_meta_p_all_cohorts"}),
    on="locus", how="left")
robust = robust[["locus"] + [c for c in robust.columns if c != "locus"]]
# S6: credible-set variants with PIP >= 0.1, plus coding variants with PIP >= 0.01 (the locus-to-gene threshold)
CODING = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant", "inframe_insertion",
          "inframe_deletion", "splice_donor_variant", "splice_acceptor_variant"}
vep = pd.read_csv(RES / "psd_finemap_vep_by_rsid.tsv", sep="\t")
vep = vep[(vep.pip_max >= 0.1) | (vep.most_severe_consequence.isin(CODING) & (vep.pip_max >= 0.01))]
# S14: MHC-free MAGMA gene sets, the top 10 per meta-analysis plus the a priori hair/epidermis/WNT sets, with robustness
# reruns and the same set in the comparator traits
gs = pd.read_csv(RES / "magma_genesets_nomhc.tsv", sep="\t")
gsr = pd.read_csv(RES / "magma_genesets_psd_robustness.tsv", sep="\t").merge(
    gs[["trait", "gene_set", "collection", "SE", "P_bonf"]], on=["trait", "gene_set"], how="left")
other = gs[~gs.trait.str.startswith("PSD")].pivot(index="gene_set", columns="trait", values="P").add_prefix("P_")
gsr = gsr.merge(other, left_on="gene_set", right_index=True, how="left").sort_values(["trait", "P"])
gsr = gsr[["trait", "gene_set", "collection", "a_priori", "NGENES", "BETA", "SE", "P", "FDR", "P_bonf", "P_noTopLocus",
           "P_noGWS", "top_gene", "top_gene_p", "n_dropped_locus", "n_dropped_gws"] + list(other.columns)]

# S1 comparator rows: the summary statistics compared with PSD (sample sizes as used in steps 02a and 04a; for BMI and
# WHRadjBMI, the largest per-variant N among HapMap3 variants)
FG = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv ; https://risteys.finngen.fi/endpoints/"
comparators = pd.DataFrame([
    ("HS_GCST90726609", "hidradenitis suppurativa", "EUR+AFR", 3941, 1435603, "GCST90726609", "Peacker et al. 2026",
     "Meta-analysis of MVP, UK Biobank and FinnGen", "https://www.ebi.ac.uk/gwas/studies/GCST90726609"),
    ("HS_FinnGen_R13", "hidradenitis suppurativa", "FIN", 1453, 474896, "FinnGen R13 L12_HIDRADENITISSUP", "Kurki et al. 2023",
     "ICD-10 L73.2, ICD-9 7058C (hospital discharge, cause-of-death registries)", FG + "L12_HIDRADENITISSUP"),
    ("ACNE_FinnGen_R13", "acne", "FIN", 4812, 474896, "FinnGen R13 L12_ACNE", "Kurki et al. 2023", "ICD-10 L70 (hospital discharge, cause-of-death registries)", FG + "L12_ACNE"),
    ("ACNE_GCST90245818", "acne", "EUR", 34422, 364991, "GCST90245818", "Teder-Laving et al. 2024",
     "Meta-analysis of three European-ancestry cohorts", "https://www.ebi.ac.uk/gwas/studies/GCST90245818"),
    ("AA_FinnGen_R13", "alopecia areata", "FIN", 994, 474896, "FinnGen R13 L12_ALOPECAREATA", "Kurki et al. 2023", "ICD-10 L63, ICD-9 7040B/7040C, ICD-8 70400 (hospital discharge, cause-of-death registries)",
     FG + "L12_ALOPECAREATA"),
    ("MPB_GCST90043616", "hair pattern 1 (no balding), men", "EUR", 66172, 140864, "GCST90043616", "Jiang et al. 2021",
     "UK Biobank data-field 2395 (men), pattern 1 vs patterns 2-4; Neale lab PHESANT coding, do not know / prefer not to answer excluded", "https://www.ebi.ac.uk/gwas/studies/GCST90043616 ; https://doi.org/10.1038/s41588-021-00954-4 (Online Methods: Neale lab pipeline for multi-categorical traits)"),
    ("MPB_Neale_2395_4", "hair pattern 4 (balding), men", "EUR", 30225, 135424, "Neale lab UK Biobank GWAS round 2, phenotype 2395_4",
     "Neale lab, round 2 (2018)", "UK Biobank data-field 2395 (men), pattern 4 vs patterns 1-3; PHESANT coding, do not know / prefer not to answer excluded", "http://www.nealelab.is/uk-biobank ; https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.male.tsv.bgz"),
    ("BMI_Pulit2019", "body mass index", "EUR", None, None, "Zenodo 1251813", "Pulit et al. 2019",
     "Quantitative, GIANT + UK Biobank, N up to 806,833", "https://zenodo.org/record/1251813"),
    ("WHRadjBMI_Pulit2019", "waist-to-hip ratio adjusted for BMI", "EUR", None, None, "Zenodo 1251813", "Pulit et al. 2019",
     "Quantitative, GIANT + UK Biobank, N up to 694,649", "https://zenodo.org/record/1251813"),
], columns=["cohort", "trait", "ancestry", "n_cases", "n_controls", "source_accession", "publication", "case_definition", "sources"])
s1 = pd.read_csv(RES / "psd_meta_cohort_qc.tsv", sep="\t").merge(pd.read_csv(BASE / "data/psd_cohort_methods.tsv", sep="\t"),
                                                                   on="cohort", how="left")
s1 = pd.concat([s1.assign(role="PSD cohort", trait="pilonidal sinus disease"), comparators.assign(role="comparator trait")],
               ignore_index=True)
s1 = s1[["role", "cohort", "trait"] + [c for c in s1.columns if c not in ("role", "cohort", "trait")]]
s1[["n_cases", "n_controls"]] = s1[["n_cases", "n_controls"]].astype("Int64")

sheets = {
    # cohort QC from the meta-analysis joined to curated per-cohort methods (data/psd_cohort_methods.tsv)
    "S1_cohorts": s1,
    "S2_lead_variants": sup,
    "S3_locus_classification": pd.read_csv(RES / "psd_meta_locus_classification.tsv", sep="\t"),
    "S4_robustness": robust,
    "S5_finemapping": pd.read_csv(RES / "psd_finemap_locus_summary.tsv", sep="\t"),
    "S6_credible_set_variants": vep,
    "S7_ldsc_h2": pd.read_csv(RES / "ldsc_v2_heritability.tsv", sep="\t"),
    "S8_ldsc_rg": pd.read_csv(RES / "ldsc_v2_genetic_correlations.tsv", sep="\t"),
    "S9_genomic_sem": pd.read_csv(RES / "genomicsem_models.tsv", sep="\t"),
    "S10_coloc_traits": pd.read_csv(RES / "coloc_psd_vs_traits.tsv", sep="\t").query("`PP.H4` >= 0.5"),
    "S11_lava": pd.read_csv(RES / "lava_bivariate_psd_noMHC.tsv", sep="\t"),
    "S12_eqtl_coloc": pd.read_csv(RES / "eqtl_coloc_psd_gtex.tsv", sep="\t").query("`PP.H4` >= 0.5"),
    "S13_locus_to_gene": l2g,
    "S14_gene_sets": gsr,
}
with pd.ExcelWriter(OUT / "supplementary_tables.xlsx") as xw:
    for name, df in sheets.items():
        df.to_excel(xw, sheet_name=name, index=False)
print(tex[:800]); print("tables written to", OUT)
