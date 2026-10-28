#!/usr/bin/env python3
"""
Figures for the PSD genetics talk (presentation/psd_genetics_talk.qmd), built from the small result tables only:
Manhattan crop (figures/psd_meta_manhattan_qq.png), LDSC rg forest plot, coloc PP.H4 heatmap, LAVA bar chart, MHC-free
MAGMA gene sets and eQTL coloc heatmap. Output: presentation/figs/*.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

BASE = Path(__file__).resolve().parent.parent
RES, OUT = BASE / "results", BASE / "presentation" / "figs"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200,
                     "savefig.bbox": "tight", "font.family": "DejaVu Sans"})
NAVY, ORANGE, GREY, RED = "#1f4e79", "#d9822b", "#9aa5b1", "#b22222"

# ---- Manhattan: top row (META_ALL Manhattan + QQ) of the existing figure ----
im = Image.open(BASE / "figures" / "psd_meta_manhattan_qq.png")
im.crop((0, 0, im.width, int(im.height * 0.47))).save(OUT / "manhattan_meta_all.png")

# ---- LDSC genetic correlations (META_EUR) ----
rg = pd.read_csv(RES / "ldsc_v2_genetic_correlations.tsv", sep="\t")
rg = rg[rg.p1 == "PSD_META_EUR"]
labels = {"HS_GCST90726609": "HS (MVP+UKB+FinnGen)", "HS_FinnGen_R13": "HS (FinnGen)",
          "ACNE_FinnGen_R13": "Acne (FinnGen)", "ACNE_GCST90245818": "Acne (Teder-Laving 2024)", "BMI_Pulit2019": "BMI",
          "WHRadjBMI_Pulit2019": "WHR adj. BMI", "BALDING4_Neale_UKB_male": "Balding, pattern 4 (UKB men)",
          "MPB_GCST90043616": "No balding, pattern 1 (UKB men)", "AA_FinnGen_R13": "Alopecia areata (FinnGen)"}
rg = rg[rg.p2.isin(labels)].assign(label=lambda d: d.p2.map(labels)).sort_values("rg")
fig, ax = plt.subplots(figsize=(9, 5.2))
y = np.arange(len(rg))
col = [NAVY if p < 0.05 / len(rg) else GREY for p in rg.p]
ax.hlines(y, rg.rg - 1.96 * rg.se, rg.rg + 1.96 * rg.se, colors=col, lw=2)
ax.scatter(rg.rg, y, c=col, s=60, zorder=3)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y, rg.label)
ax.set_xlabel("Genetic correlation with PSD (LDSC rg, 95% CI)")
ax.set_xlim(-1, 1)
for yi, (r, p) in enumerate(zip(rg.rg, rg.p)):
    ax.text(1.02, yi, f"{r:+.2f}  (p={p:.0e})", va="center", ha="left", fontsize=11, transform=ax.get_yaxis_transform())
ax.set_title("Blue: Bonferroni-significant", fontsize=11, loc="left", color=NAVY)
fig.savefig(OUT / "ldsc_rg.png"); plt.close(fig)

# ---- coloc heatmap (META_EUR, non-MHC PSD loci) ----
co = pd.read_csv(RES / "coloc_psd_vs_traits.tsv", sep="\t")
co = co[(co.psd_set == "PSD_META_EUR") & (co.region_source == "PSD") & (~co.is_mhc)]
traits = ["HS", "HS_FG", "ACNE_FG", "ACNE_TL", "BALDP1", "BMI", "WHRadjBMI"]
tlab = ["HS", "HS\n(FinnGen)", "Acne\n(FinnGen)", "Acne\n(T-L 2024)", "Balding", "BMI", "WHR\nadj BMI"]
co["locus"] = co.label.where(~co.label.str.startswith("chr"), co.label.str.replace(r"^chr(\d+):(\d+)$",
              lambda m: f"chr{m.group(1)}:{int(m.group(2)) / 1e6:.1f}Mb", regex=True))
co["locus"] = co.locus + np.where(co.psd_status == "KNOWN", " *", "")
h4 = co.pivot_table(index="locus", columns="trait", values="PP.H4").reindex(columns=traits)
h3 = co.pivot_table(index="locus", columns="trait", values="PP.H3").reindex(columns=traits)
order = (h4 > 0.8).sum(axis=1).sort_values(ascending=False).index
order = sorted(order, key=lambda l: (-(h4.loc[l] > 0.8).sum(), -(h4.loc[l, ["HS", "HS_FG"]] > 0.8).sum(), l))
h4, h3 = h4.loc[order], h3.loc[order]
fig, ax = plt.subplots(figsize=(8.2, 9))
m = ax.imshow(h4.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i in range(h4.shape[0]):
    for j in range(h4.shape[1]):
        v4, v3 = h4.iat[i, j], h3.iat[i, j]
        if v4 > 0.8:
            ax.text(j, i, f"{v4:.2f}", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        elif v3 > 0.8:
            ax.text(j, i, "×", ha="center", va="center", color=RED, fontsize=13)
ax.set_xticks(range(len(traits)), tlab, fontsize=10)
ax.set_yticks(range(len(order)), order, fontsize=10)
ax.xaxis.tick_top()
cb = fig.colorbar(m, ax=ax, fraction=0.04, pad=0.02); cb.set_label("PP.H4 (shared causal variant)")
ax.set_xlabel("number = PP.H4 > 0.8     × = distinct causal variants (PP.H3 > 0.8)     * = previously reported locus", fontsize=9.5)
fig.savefig(OUT / "coloc_heatmap.png"); plt.close(fig)

# ---- LAVA local rg (non-MHC, BH FDR < 0.05) ----
lv = pd.read_csv(RES / "lava_bivariate_psd_noMHC.tsv", sep="\t")
lmap = {"HS": "HS", "HS_FG": "HS (FinnGen)", "ACNE_FG": "Acne (FinnGen)", "ACNE_TL": "Acne (T-L 2024)",
        "BALDP1": "No balding (pattern 1)", "BMI": "BMI", "WHRadjBMI": "WHR adj. BMI"}
lv = lv[lv.phen1.isin(lmap)]
s = lv.assign(sig=lv.fdr < 0.05, pos=lv.rho > 0).groupby("phen1").agg(
    tested=("p", "size"), pos=("pos", lambda x: (x & lv.loc[x.index, "fdr"].lt(0.05)).sum()),
    neg=("pos", lambda x: ((~x) & lv.loc[x.index, "fdr"].lt(0.05)).sum())).reindex(list(lmap)).dropna()
fig, ax = plt.subplots(figsize=(8.5, 4.6))
x = np.arange(len(s))
ax.bar(x, s.pos, color=NAVY, label="positive local rg")
ax.bar(x, s.neg, bottom=s.pos, color=ORANGE, label="negative local rg")
for xi, (p, n, t) in enumerate(zip(s.pos, s.neg, s.tested)):
    ax.text(xi, p + n + 0.4, f"{int(p + n)}/{int(t)}", ha="center", fontsize=11)
ax.set_xticks(x, [lmap[k] for k in s.index], rotation=25, ha="right")
ax.set_ylabel("Blocks with FDR < 0.05")
ax.legend(frameon=False, loc="upper right")
fig.savefig(OUT / "lava_counts.png"); plt.close(fig)

# ---- MHC-free MAGMA gene sets: top PSD_META_ALL sets, with META_EUR and the comparator traits ----
gs = pd.read_csv(RES / "magma_genesets_nomhc.tsv", sep="\t")
top = gs[gs.trait == "PSD_META_ALL"].nsmallest(8, "P").gene_set.tolist()
lp = gs[gs.gene_set.isin(top)].assign(mlp=lambda d: -np.log10(d.P)).pivot(index="gene_set", columns="trait", values="mlp").reindex(top[::-1])
fig, ax = plt.subplots(figsize=(10, 4.8))
y = np.arange(len(top))
ax.barh(y + 0.2, lp.PSD_META_ALL, height=0.4, color=NAVY, label="PSD, META_ALL")
ax.barh(y - 0.2, lp.PSD_META_EUR, height=0.4, color=GREY, label="PSD, META_EUR")
ax.axvline(-np.log10(0.05 / gs[gs.trait == "PSD_META_ALL"].shape[0]), color=RED, ls="--", lw=1.2, label="Bonferroni (6,717 sets)")
ax.set_yticks(y, [("GO: " if s.startswith("GOBP_") else "Reactome: ") + s.split("_", 1)[1].replace("_", " ").lower()[:48] for s in lp.index], fontsize=10.5)
ax.set_xlabel("−log10 p (MAGMA competitive test, MHC genes removed)")
ax.legend(frameon=False, loc="lower right", fontsize=10)
fig.savefig(OUT / "genesets.png"); plt.close(fig)

# ---- GTEx eQTL coloc (META_EUR): locus -> gene pairs with PP.H4 >= 0.6 in any tissue ----
eq = pd.read_csv(RES / "eqtl_coloc_psd_gtex.tsv", sep="\t")
eq = eq[eq.psd_set == "PSD_META_EUR"].copy()
# Ensembl GRCh38 (checked 2026-10-21): only TANAR and RMDN2-AS1 have official symbols; the others are unnamed novel transcripts
lnc = {"ENSG00000287727": "unnamed antisense lncRNA", "ENSG00000223991": "TANAR (lncRNA)", "ENSG00000236838": "unnamed antisense lncRNA",
       "ENSG00000232821": "lncRNA overlapping TWIST1", "ENSG00000235848": "RMDN2-AS1 (lncRNA)"}
eq["gene"] = eq.gene_id.map(lnc).fillna(eq.gene_name)
eq["locus"] = eq.label.str.replace(r"^chr(\d+):(\d+)$", lambda m: f"chr{m.group(1)}:{int(m.group(2)) / 1e6:.1f}Mb", regex=True)
eq["pair"] = eq.locus + "  →  " + eq.gene
tis = ["Skin_sun_exposed", "Skin_not_sun_exposed", "Fibroblast_cultured", "Adipose_subcutaneous", "Adipose_visceral"]
tlab = ["Skin\n(sun-exp.)", "Skin\n(not sun-exp.)", "Fibroblasts", "Adipose\n(subcut.)", "Adipose\n(visceral)"]
keep = eq.groupby("pair")["PP.H4"].max().loc[lambda s: s >= 0.6].sort_values(ascending=False).index
h4 = eq[eq.pair.isin(keep)].pivot_table(index="pair", columns="tissue", values="PP.H4").reindex(index=keep, columns=tis)
h3 = eq[eq.pair.isin(keep)].pivot_table(index="pair", columns="tissue", values="PP.H3").reindex(index=keep, columns=tis)
fig, ax = plt.subplots(figsize=(9, 4.6))
m = ax.imshow(h4.fillna(0).values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i in range(h4.shape[0]):
    for j in range(h4.shape[1]):
        v4, v3 = h4.iat[i, j], h3.iat[i, j]
        if pd.isna(v4):
            ax.text(j, i, "·", ha="center", va="center", color=GREY, fontsize=12)
        elif v4 >= 0.5:
            ax.text(j, i, f"{v4:.2f}", ha="center", va="center", color="white" if v4 >= 0.7 else "black", fontsize=10, fontweight="bold")
        elif v3 >= 0.8:
            ax.text(j, i, "×", ha="center", va="center", color=RED, fontsize=14)
ax.set_xticks(range(len(tis)), tlab, fontsize=10)
ax.set_yticks(range(len(keep)), keep, fontsize=11)
ax.xaxis.tick_top()
cb = fig.colorbar(m, ax=ax, fraction=0.035, pad=0.02); cb.set_label("PP.H4", fontsize=10)
ax.set_xlabel("× = distinct causal variants (PP.H3 > 0.8)     · = gene not tested (eQTL p > 1e-4)", fontsize=9.5)
fig.savefig(OUT / "eqtl_coloc.png"); plt.close(fig)
print("figures written to", OUT)
