#!/usr/bin/env python3
"""Figures for the manuscript ($PSD_PAPER/manuscript.tex), built locally from the small result tables in results/.
Figure 1 (Miami plot) and Supplementary Figure S1 (QQ) come from scripts/11a_paper_manhattan_qq.py (cluster) and 11b_paper_miami.py.
Outputs: $PSD_PAPER/figures/fig2_shared_genetics.png, fig3_coloc.png, fig4_effector_genes.png, figS2_cohort_forest.png
"""
import os, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = Path(__file__).resolve().parents[1]
PAPER = Path(os.environ.get("PSD_PAPER", BASE / "paper"))  # manuscript directory (kept off the public repo)
RES, OUT = Path(os.environ.get("PSD_RESULTS", BASE / "results")), PAPER / "figures"

def save(fig, name):
    """PNG for review drafts, PDF (vector, editable text) for submission."""
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "pdf.fonttype": 42, "savefig.bbox": "tight", "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.titlelocation": "left"})
NAVY, ORANGE, GREY, RED = "#1f4e79", "#d9822b", "#9aa5b1", "#b22222"

def locus_name(label):
    m = re.match(r"^chr(\d+):(\d+)$", str(label))
    return f"chr{m.group(1)}:{int(m.group(2)) / 1e6:.1f} Mb" if m else label

# ---------------- Figure 2: genome-wide and local genetic correlation ----------------
rg = pd.read_csv(RES / "ldsc_v2_genetic_correlations.tsv", sep="\t")
rg = rg[rg.p1 == "PSD_META_EUR"]
lab = {"HS_GCST90726609": "HS", "HS_FinnGen_R13": "HS (FinnGen)", "ACNE_FinnGen_R13": "Acne (FinnGen)",
       "ACNE_GCST90245818": "Acne (Teder-Laving)", "MPB_GCST90043616": "No balding (pattern 1)",
       "BALDING4_Neale_UKB_male": "Balding (pattern 4)", "BMI_Pulit2019": "BMI", "WHRadjBMI_Pulit2019": "WHR adj. BMI",
       "AA_FinnGen_R13": "Alopecia areata"}
rg = rg[rg.p2.isin(lab)].assign(label=lambda d: d.p2.map(lab)).sort_values("rg")
sem = pd.read_csv(RES / "genomicsem_models.tsv", sep="\t")
sem = sem[(sem.op == "~") & (sem.lhs == "PSD")]
slab = {"HS": "HS", "BMI": "BMI", "ACNE_FG": "Acne (FinnGen)", "WHRadjBMI": "WHR adj. BMI", "BALDP1": "No balding (pattern 1)"}
lv = pd.read_csv(RES / "lava_bivariate_psd_noMHC.tsv", sep="\t")
lmap = {"HS": "HS", "ACNE_FG": "Acne (FinnGen)", "ACNE_TL": "Acne (Teder-Laving)", "BALDP1": "No balding (pattern 1)",
        "BMI": "BMI", "WHRadjBMI": "WHR adj. BMI"}
lv = lv[lv.phen1.isin(lmap)].assign(sig=lambda d: d.fdr < 0.05)
cnt = lv.groupby("phen1").apply(lambda x: pd.Series({"tested": len(x), "pos": (x.sig & (x.rho > 0)).sum(),
                                                       "neg": (x.sig & (x.rho < 0)).sum()})).reindex(list(lmap))

fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 1.25})
ax = axes[0]; y = np.arange(len(rg)); col = [NAVY if p < 0.05 / len(rg) else GREY for p in rg.p]
ax.hlines(y, rg.rg - 1.96 * rg.se, rg.rg + 1.96 * rg.se, colors=col, lw=1.4)
ax.scatter(rg.rg, y, c=col, s=14, zorder=3); ax.axvline(0, color="black", lw=0.6)
ax.set_yticks(y, rg.label); ax.set_xlim(-1, 1.1); ax.set_xlabel("Genetic correlation (rg, 95% CI)")
ax.set_title("(a)  Genome-wide\n      (LDSC)")
ax = axes[1]
rows = [("M2a_PSD_on_HS_BMI", "HS"), ("M2a_PSD_on_HS_BMI", "BMI"), ("M2b_PSD_on_all", "HS"), ("M2b_PSD_on_all", "ACNE_FG"),
        ("M2b_PSD_on_all", "BMI"), ("M2b_PSD_on_all", "WHRadjBMI"), ("M2b_PSD_on_all", "BALDP1")]
yy = [6, 5, 3.5, 2.5, 1.5, 0.5, -0.5]
for (m, t), yi in zip(rows, yy):
    r = sem[(sem.model == m) & (sem.rhs == t)].iloc[0]
    c = NAVY if r.p_value < 0.05 else GREY
    ax.hlines(yi, r.STD_Genotype - 1.96 * r.STD_Genotype_SE, r.STD_Genotype + 1.96 * r.STD_Genotype_SE, colors=c, lw=1.4)
    ax.scatter(r.STD_Genotype, yi, color=c, s=14, zorder=3)
ax.set_yticks(yy, [slab[t] for _, t in rows]); ax.axvline(0, color="black", lw=0.6)
ax.axhline(4.4, color=GREY, lw=0.5, ls=":")
ax.text(0.98, 6.6, "Model 1", transform=ax.get_yaxis_transform(), fontsize=6.5, color="#555", va="bottom", ha="right")
ax.text(0.98, 4.1, "Model 2", transform=ax.get_yaxis_transform(), fontsize=6.5, color="#555", va="bottom", ha="right")
ax.set_ylim(-1, 7.2)
ax.set_xlabel("Standardised β on PSD (95% CI)"); ax.set_title("(b)  Joint model\n      (Genomic SEM)")
ax = axes[2]; x = np.arange(len(cnt))
ax.barh(x, cnt.pos, color=NAVY, label="positive"); ax.barh(x, cnt.neg, left=cnt.pos, color=ORANGE, label="negative")
for xi, (p, n, t) in enumerate(zip(cnt.pos, cnt.neg, cnt.tested)):
    ax.text(p + n + 0.4, xi, f"{int(p + n)}/{int(t)}", va="center", fontsize=7)
ax.set_yticks(x, [lmap[k] for k in cnt.index]); ax.invert_yaxis(); ax.set_xlim(0, 26)
ax.set_xlabel("Blocks with local rg (FDR < 0.05)"); ax.set_title("(c)  Local\n      (LAVA)")
ax.legend(frameon=False, fontsize=6.5, loc="upper right", bbox_to_anchor=(1.05, 0.62))
save(fig, "fig2_shared_genetics")

# ---------------- Figure 3: shared causal variants (coloc) ----------------
co = pd.read_csv(RES / "coloc_psd_vs_traits.tsv", sep="\t")
co = co[(co.psd_set == "PSD_META_EUR") & (co.region_source == "PSD") & (~co.is_mhc.astype(bool))]
traits = ["HS", "HS_FG", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI"]
tlab = ["HS", "HS\n(FinnGen)", "Acne\n(FinnGen)", "Acne\n(Teder-\nLaving)", "Balding", "BMI", "WHR\nadj. BMI"]
co["locus"] = co.label.map(locus_name) + np.where(co.psd_status == "KNOWN", " *", "")
h4 = co.pivot_table(index="locus", columns="trait", values="PP.H4").reindex(columns=traits)
h3 = co.pivot_table(index="locus", columns="trait", values="PP.H3").reindex(columns=traits)
order = sorted(h4.index, key=lambda l: (-int((h4.loc[l, ["HS", "HS_FG", "ACNE_FG", "ACNE_TL"]] > 0.8).any()),
                                        -int(h4.loc[l, "BALDP1"] > 0.8), -int((h4.loc[l] > 0.8).sum()), l))
h4, h3 = h4.loc[order], h3.loc[order]
fig, ax = plt.subplots(figsize=(4.6, 6.4))
m = ax.imshow(h4.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i in range(h4.shape[0]):
    for j in range(h4.shape[1]):
        v4, v3 = h4.iat[i, j], h3.iat[i, j]
        if v4 > 0.8:
            ax.text(j, i, f"{v4:.2f}", ha="center", va="center", color="white", fontsize=6.5, fontweight="bold")
        elif v3 > 0.8:
            ax.text(j, i, "×", ha="center", va="center", color=RED, fontsize=9)
ax.set_xticks(range(len(traits)), tlab, fontsize=7); ax.xaxis.tick_top()
ax.set_yticks(range(len(order)), order, fontsize=7)
for t in ax.get_yticklabels():
    if not t.get_text().startswith("chr"):
        t.set_fontstyle("italic")
cb = fig.colorbar(m, ax=ax, fraction=0.05, pad=0.03); cb.set_label("PP.H4 (shared causal variant)")
ax.set_xlabel("number: PP.H4 > 0.8   ×: distinct causal variants (PP.H3 > 0.8)\n*: previously reported PSD locus", fontsize=7)
save(fig, "fig3_coloc")

# ---------------- Figure 4: candidate effector genes ----------------
eq = pd.read_csv(RES / "eqtl_coloc_psd_gtex.tsv", sep="\t")
eq = eq[eq.psd_set == "PSD_META_EUR"].copy()
lnc = {"ENSG00000287727": "lncRNA, antisense", "ENSG00000223991": "TANAR (lncRNA)", "ENSG00000236838": "lncRNA, antisense",
       "ENSG00000232821": "lncRNA (TWIST1)", "ENSG00000235848": "RMDN2-AS1 (lncRNA)"}
eq["gene"] = eq.gene_id.map(lnc).fillna(eq.gene_name)
eq["pair"] = eq.label.map(locus_name) + " → " + eq.gene
tis = ["Skin_sun_exposed", "Skin_not_sun_exposed", "Fibroblast_cultured", "Adipose_subcutaneous", "Adipose_visceral"]
tl = ["Skin\n(sun-exp.)", "Skin\n(not sun-exp.)", "Fibroblasts", "Adipose\n(subcut.)", "Adipose\n(visceral)"]
keep = eq.groupby("pair")["PP.H4"].max().loc[lambda s: s >= 0.6].sort_values(ascending=False).index
e4 = eq[eq.pair.isin(keep)].pivot_table(index="pair", columns="tissue", values="PP.H4").reindex(index=keep, columns=tis)
e3 = eq[eq.pair.isin(keep)].pivot_table(index="pair", columns="tissue", values="PP.H3").reindex(index=keep, columns=tis)
direction = eq[eq.pair.isin(keep)].sort_values("PP.H4").groupby("pair").risk_allele_expression.last()

l2g = pd.read_csv(RES / "psd_locus_to_gene.tsv", sep="\t")
lines = ["nearest", "MAGMA", "coding", "eQTL", "L2G"]
def evid(row):
    g0 = str(row.prioritised_gene).split("/")[0].strip()
    ent = {e.split(":")[0].strip(): e.split(":", 1)[1] for e in str(row.evidence).split("; ") if ":" in e}
    txt = ent.get(g0, "")
    vals = {"nearest": 1 if "nearest" in txt else 0, "MAGMA": 1 if "MAGMA" in txt else 0,
            "coding": 2 if re.search(r"coding PIP (0\.[1-9]|1)", txt) else (1 if "coding" in txt else 0),
            "eQTL": 2 if re.search(r"eQTL coloc (0\.[89]|1)", txt) else (1 if "eQTL" in txt else 0), "L2G": 1 if "L2G" in txt else 0}
    return pd.Series(vals)
ev = l2g.apply(evid, axis=1)
def short(g, near):
    g = str(g)
    if "antisense" in g or g.startswith("ENSG"):
        return f"lncRNA ({near} locus)"
    return g.replace(" (lncRNA)", " (lncRNA)").split("/")[0] if "lncRNA" in g else g
ev.index = [short(g, n) for g, n in zip(l2g.prioritised_gene, l2g.nearest_gene)]
tier_order = {"functional (protein-coding)": 0, "functional (lncRNA only)": 1, "supportive (2+ lines)": 2, "single line": 3}
o = np.lexsort((-l2g.evidence_score.values, l2g.locus_evidence_tier.map(tier_order).values))
ev = ev.iloc[o]

fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.6), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 1.1})
ax = axes[0]
ax.imshow(e4.fillna(0).values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i in range(e4.shape[0]):
    for j in range(e4.shape[1]):
        v4, v3 = e4.iat[i, j], e3.iat[i, j]
        if pd.isna(v4):
            ax.text(j, i, "·", ha="center", va="center", color=GREY, fontsize=8)
        elif v4 >= 0.5:
            ax.text(j, i, f"{v4:.2f}", ha="center", va="center", color="white" if v4 >= 0.7 else "black", fontsize=6.5, fontweight="bold")
        elif v3 >= 0.8:
            ax.text(j, i, "×", ha="center", va="center", color=RED, fontsize=9)
ax.set_xticks(range(len(tis)), [t.replace("\n", " ") for t in tl], fontsize=6.5, rotation=45, ha="left"); ax.xaxis.tick_top()
ax.set_yticks(range(len(keep)), [f"{p} ({'↑' if direction[p] == 'up' else '↓'})" for p in keep], fontsize=7)
ax.set_xlabel("PP.H4 shown if ≥ 0.5; ×: distinct variants; ·: not tested\n↑/↓: risk allele raises/lowers expression", fontsize=6.5)
ax.set_title("(a)  eQTL colocalisation (GTEx v8)", pad=48)
ax = axes[1]
cmap = matplotlib.colors.ListedColormap(["#f2f2f2", "#9ecae1", NAVY])
ax.imshow(ev.values, cmap=cmap, vmin=0, vmax=2, aspect="auto")
ax.set_xticks(range(len(lines)), lines, fontsize=6.5, rotation=45, ha="left"); ax.xaxis.tick_top()
ax.set_yticks(range(len(ev)), ev.index, fontsize=6.5)
for t in ax.get_yticklabels():
    t.set_fontstyle("italic")
bounds = np.cumsum(l2g.locus_evidence_tier.map(tier_order).value_counts().sort_index().values)[:-1]
for b in bounds:
    ax.axhline(b - 0.5, color="black", lw=0.8)
ax.set_title("(b)  Locus-to-gene evidence", pad=48)
ax.legend([Line2D([], [], marker="s", ls="", color=c, ms=6) for c in ("#9ecae1", NAVY)], ["weak", "strong"],
          frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.02, 1))
save(fig, "fig4_effector_genes")

# ---------------- Supplementary Figure S2: per-cohort effects at the 25 lead variants ----------------
sup = pd.read_csv(RES / "psd_lead_variants_supplementary.tsv", sep="\t")
sup = sup[sup.meta == "META_EUR"].reset_index(drop=True)
def parse(s):
    m = re.match(r"([\d.]+) \(([\d.]+)-([\d.]+)\)", str(s)); return (float(m[1]), float(m[2]), float(m[3])) if m else (np.nan,) * 3
cols = {"MVP_EUR": ("MVP", "#2b8cbe"), "UKB_EUR": ("UKB", "#7bccc4"), "FinnGen_R13": ("FinnGen", "#e6550d"), "fixed": ("Meta (fixed)", "black")}
fig, ax = plt.subplots(figsize=(4.6, 7.2))
for k, (c, (name, color)) in enumerate(cols.items()):
    vals = sup["or_fixed_95ci" if c == "fixed" else f"{c}_or_95ci"].map(parse)
    for i, (o, lo, hi) in enumerate(vals):
        if np.isnan(o): continue
        yv = i + (k - 1.5) * 0.18
        ax.hlines(yv, lo, hi, color=color, lw=0.9); ax.scatter(o, yv, color=color, s=9 if c != "fixed" else 16,
                                                              marker="o" if c != "fixed" else "D", zorder=3)
ax.set_xscale("log"); ax.axvline(1, color="black", lw=0.6)
ax.set_xticks([0.4, 0.6, 0.8, 1, 1.25, 1.5], ["0.4", "0.6", "0.8", "1", "1.25", "1.5"]); ax.minorticks_off()
ax.set_yticks(range(len(sup)), [f"{n} ({r})" for n, r in zip(sup.nearest_gene.fillna(sup.locus), sup.rsid)], fontsize=6.5)
ax.invert_yaxis(); ax.set_xlabel("Odds ratio (95% CI), effect allele as in Table 1")
ax.legend([Line2D([], [], color=c, marker="D" if n.startswith("Meta") else "o", ms=4, lw=0.9) for n, c in cols.values()],
          [n for n, _ in cols.values()], frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4)
save(fig, "figS2_cohort_forest")
print("figures written to", OUT)
