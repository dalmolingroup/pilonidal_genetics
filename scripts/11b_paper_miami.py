#!/usr/bin/env python3
"""Paper Figure 1: Miami plot of META_ALL (top) and META_EUR (bottom) from $PSD_PAPER/data/miami_thinned.tsv.gz
(written on the cluster by scripts/11a_paper_manhattan_qq.py: variants in >= 2 cohorts, p > 0.01 thinned to 10%). GRCh38.
The 25 non-MHC loci are labelled with their prioritised gene (red: novel), spread apart with leader lines.
Output: $PSD_PAPER/figures/fig1_miami.{png,pdf}
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
PAPER = Path(os.environ.get("PSD_PAPER", BASE / "paper"))  # manuscript directory (kept off the public repo)
(PAPER / "figures").mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
GWS, MHC, YMAX = 5e-8, (25_000_000, 34_000_000), 27

d = pd.read_csv(PAPER / "data/miami_thinned.tsv.gz", sep="\t")
l2g = pd.read_csv(Path(os.environ.get("PSD_RESULTS", BASE / "results")) / "psd_locus_to_gene.tsv", sep="\t")
loc = l2g.locus.str.extract(r"chr(\d+):(\d+)").astype(int)
def label(p, n):
    p = str(p)
    if "antisense" in p or p.startswith("ENSG"):
        return n
    return p.split("/")[0].replace(" (lncRNA)", "")
lab = pd.DataFrame({"chr": loc[0], "pos": loc[1], "gene": [label(p, n) for p, n in zip(l2g.prioritised_gene, l2g.nearest_gene)],
                    "novel": l2g.status.eq("NOVEL")})

sizes = d.groupby("chr").pos.max()
offset = pd.Series(np.concatenate([[0], np.cumsum(sizes.values[:-1] + 8e6)]), index=sizes.index)
total = offset.iloc[-1] + sizes.iloc[-1]

fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True, gridspec_kw={"hspace": 0.04, "height_ratios": [1.25, 1]})
for ax, meta, sign in zip(axes, ("ALL", "EUR"), (1, -1)):
    df = d[d.meta == meta]
    x = df.pos + df.chr.map(offset); y = -np.log10(df.p.clip(lower=1e-300)) * sign
    col = np.where(df.chr % 2 == 0, "#9aa7b8", "#2b4a6f")
    col = np.where((df.chr == 6) & df.pos.between(*MHC), "#d9d9d9", col)
    ax.scatter(x, y, c=col, s=1.2, lw=0, rasterized=True)
    ax.axhline(-np.log10(GWS) * sign, color="#c0392b", lw=0.6, ls="--")
    ax.set_ylabel(f"META_{meta}\n−log₁₀ p", fontsize=8)
axes[0].set_ylim(0, YMAX + 9); axes[0].set_yticks(range(0, 26, 5))
axes[1].set_ylim(-YMAX, 0); axes[1].set_yticks([-25, -20, -15, -10, -5, 0], ["25", "20", "15", "10", "5", "0"])
centers = offset + sizes / 2
axes[1].set_xticks(centers.values, [str(c) if c < 17 or c % 2 else "" for c in centers.index], fontsize=7)
axes[1].set_xlabel("Chromosome (GRCh38)")
axes[1].set_xlim(-1e7, total + 1e7)

# labels: one row above the tallest peak band, x spread so neighbours are at least `gap` apart, leader lines to the peaks
top = d[d.meta == "ALL"]
lab["x"] = lab.pos + lab.chr.map(offset)
lab["peak"] = [(-np.log10(w.p.min()) if len(w) else 8) for w in
               (top[(top.chr == r.chr) & top.pos.between(r.pos - 500_000, r.pos + 500_000)] for r in lab.itertuples())]
lab = lab.sort_values("x").reset_index(drop=True)
gap = total * 0.026
xs = lab.x.values.astype(float).copy()
for _ in range(200):                                   # simple 1-D repulsion
    moved = False
    for i in range(1, len(xs)):
        if xs[i] - xs[i - 1] < gap:
            shift = (gap - (xs[i] - xs[i - 1])) / 2; xs[i - 1] -= shift; xs[i] += shift; moved = True
    if not moved: break
ylab = YMAX + 0.5
for r, xl in zip(lab.itertuples(), xs):
    c = "#c0392b" if r.novel else "black"
    ax = axes[0]
    ax.plot([r.x, r.x, xl], [r.peak + 0.4, ylab - 1.2, ylab - 0.3], color="#999", lw=0.4)
    ax.text(xl, ylab, r.gene, rotation=90, ha="center", va="bottom", fontsize=6, fontstyle="italic", color=c)
for ext in ("png", "pdf"):                             # scatter points are rasterised, text stays editable
    fig.savefig(PAPER / f"figures/fig1_miami.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("fig1 written")
