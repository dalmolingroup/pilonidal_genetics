#!/usr/bin/env python3
"""Paper Figure 1b-c (cluster): Miami plot of META_ALL (top) and META_EUR (bottom), with the 25 non-MHC loci labelled,
and QQ plots for the supplement. Reads data/processed/meta/psd_meta_{ALL,EUR}.tsv.gz (variants in >= 2 cohorts only, as for
discovery); variants with p > 0.01 are thinned to 10% for plotting. GRCh38.
Outputs: $PSD_PAPER/data/miami_thinned.tsv.gz (re-plotted locally by scripts/11b_paper_miami.py), $PSD_PAPER/figures/figS1_qq.png
"""
import os
from pathlib import Path
import duckdb, numpy as np, pandas as pd
from scipy.stats import chi2
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
PAPER = Path(os.environ.get("PSD_PAPER", BASE / "paper"))  # manuscript directory (kept off the public repo)
OUT = PAPER / "figures"; OUT.mkdir(parents=True, exist_ok=True); (PAPER / "data").mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
GWS, MHC = 5e-8, (25_000_000, 34_000_000)
con = duckdb.connect(); con.execute(f"SET threads={os.environ.get('NCPU', 4)}")

def load(meta):
    f = BASE / f"data/processed/meta/psd_meta_{meta}.tsv.gz"
    return con.execute(f"""SELECT CAST(chr AS INTEGER) AS chr, pos, p_meta AS p FROM read_csv_auto('{f}', delim='\t', header=true)
                           WHERE k >= 2 AND TRY_CAST(chr AS INTEGER) BETWEEN 1 AND 22
                             AND (p_meta < 0.01 OR hash(chr || ':' || pos) % 10 = 0)""").df()

l2g = pd.read_csv(BASE / "results/psd_locus_to_gene.tsv", sep="\t")
loc = l2g.locus.str.extract(r"chr(\d+):(\d+)").astype(int)
labels = pd.DataFrame({"chr": loc[0], "pos": loc[1],
                       "gene": [p.split("/")[0].replace(" (lncRNA)", "") if not p.startswith("ENSG") else n
                                for p, n in zip(l2g.prioritised_gene, l2g.nearest_gene)],
                       "novel": l2g.status.eq("NOVEL")})

d = {m: load(m) for m in ("ALL", "EUR")}
pd.concat([v.assign(meta=m) for m, v in d.items()]).to_csv(OUT.parent / "data/miami_thinned.tsv.gz", sep="\t", index=False)
sizes = pd.concat(d.values()).groupby("chr").pos.max()
offset = pd.Series(np.concatenate([[0], np.cumsum(sizes.values[:-1] + 5e6)]), index=sizes.index)
centers = offset + sizes / 2

fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True, gridspec_kw={"hspace": 0.05})
for ax, (meta, df), sign in zip(axes, d.items(), (1, -1)):
    x = df.pos + df.chr.map(offset); y = -np.log10(df.p.clip(lower=1e-300)) * sign
    col = np.where(df.chr % 2 == 0, "#9aa7b8", "#2b4a6f")
    mhc = (df.chr == 6) & df.pos.between(*MHC)
    col = np.where(mhc, "#d0d0d0", col)
    ax.scatter(x, y, c=col, s=1.2, lw=0, rasterized=True)
    ax.axhline(-np.log10(GWS) * sign, color="#c0392b", lw=0.6, ls="--")
    ax.set_ylabel(f"META_{meta}\n−log₁₀ p")
    ax.set_ylim((0, 27) if sign > 0 else (-27, 0))
    if sign < 0:
        ax.set_yticks([-25, -20, -15, -10, -5, 0], ["25", "20", "15", "10", "5", "0"])
axes[1].set_xticks(centers.values, [str(c) if c < 19 or c % 2 else "" for c in centers.index], fontsize=7)
axes[1].set_xlabel("Chromosome (GRCh38)")
top = d["ALL"]
for r in labels.itertuples():
    w = top[(top.chr == r.chr) & (top.pos.between(r.pos - 500_000, r.pos + 500_000))]
    yv = -np.log10(w.p.min()) if len(w) else 8
    axes[0].annotate(r.gene, (r.pos + offset[r.chr], yv), xytext=(0, 3), textcoords="offset points", rotation=90,
                     ha="center", va="bottom", fontsize=6, fontstyle="italic", color="#c0392b" if r.novel else "black")
fig.savefig(OUT / "fig1_miami.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# QQ plots (supplement), from the full distribution (no thinning)
fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0))
for ax, meta in zip(axes, ("ALL", "EUR")):
    f = BASE / f"data/processed/meta/psd_meta_{meta}.tsv.gz"
    p = np.sort(con.execute(f"SELECT p_meta FROM read_csv_auto('{f}', delim='\t', header=true) WHERE k >= 2").df().p_meta.values)
    n = len(p); exp = -np.log10((np.arange(1, n + 1) - 0.5) / n); obs = -np.log10(p.clip(1e-300))
    keep = (obs > 2) | (np.arange(n) % 200 == 0)
    lam = chi2.isf(np.median(p), 1) / chi2.isf(0.5, 1)
    ax.scatter(exp[keep], obs[keep], s=2, c="#2b4a6f", lw=0, rasterized=True)
    ax.plot([0, exp.max()], [0, exp.max()], color="#c0392b", lw=0.6)
    ax.set_title(f"META_{meta} (λGC = {lam:.3f}; {n:,} variants)", fontsize=8)
    ax.set_xlabel("Expected −log₁₀ p"); ax.set_ylabel("Observed −log₁₀ p")
fig.tight_layout(); fig.savefig(OUT / "figS1_qq.png", dpi=300); plt.close(fig)
print("MANHATTAN_DONE")
