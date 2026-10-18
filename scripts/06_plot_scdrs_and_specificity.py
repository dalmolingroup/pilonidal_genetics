#!/usr/bin/env python3
"""
Step 6: Summarize and visualize scDRS cell-type enrichment (association & heterogeneity p-values,
mean normalized scDRS score, and proportion of disease-relevant cells) alongside single-gene
expression specificity and GTEx eQTL regulation across skin, hair-follicle, and mesenchymal compartments.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)


def main():
    stat_files = sorted(RESULTS_DIR.glob("scdrs_group_stats_*.tsv"))
    if not stat_files:
        print("No scdrs_group_stats_*.tsv files found yet.")
        return

    dfs = [pd.read_csv(f, sep="\t") for f in stat_files]
    df_all = pd.concat(dfs, ignore_index=True)
    df_all = df_all[~df_all["trait"].astype(str).str.startswith("PSD_Module_")].copy()
    df_all["neg_log10_assoc_p"] = -np.log10(df_all["assoc_mcp"].clip(lower=1e-4))
    df_all["neg_log10_hetero_p"] = -np.log10(df_all["hetero_mcp"].clip(lower=1e-4))
    out_combined = RESULTS_DIR / "scdrs_all_atlases_combined_summary.tsv"
    df_all.sort_values(["assoc_mcp", "mean_zscore"], ascending=[True, False]).to_csv(out_combined, sep="\t", index=False)
    print(f"Saved combined scDRS summary ({len(df_all)} rows) -> {out_combined}")

    # Print top enriched cell types per trait and atlas
    top_hits = df_all[df_all["assoc_mcp"] < 0.10].sort_values(["assoc_mcp", "mean_zscore"], ascending=[True, False])
    print(f"\n=== Top scDRS Enriched Cell Types / Subclusters (assoc_mcp < 0.10: {len(top_hits)} hits) ===")
    cols_show = ["atlas", "trait", "group_col", "group_level", "n_cell", "assoc_mcp", "hetero_mcp", "mean_zscore", "prop_p_lt_005", "n_fdr_0.1"]
    print(top_hits[cols_show].head(40).to_string(index=False))

    # Plot 1: Heatmap of scDRS mean_zscore and significance across cell types for key atlases
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    for atlas_name in df_all["atlas"].unique():
        sub = df_all[df_all["atlas"] == atlas_name].copy()
        primary_gcol = sub["group_col"].iloc[0]
        sub = sub[sub["group_col"] == primary_gcol]
        if sub.empty:
            continue
        piv_z = sub.pivot(index="group_level", columns="trait", values="assoc_mcz").fillna(0.0)
        piv_p = sub.pivot(index="group_level", columns="trait", values="assoc_mcp").fillna(1.0)
        # Annotate with asterisks for assoc_mcp < 0.05 (*) and < 0.01 (**)
        annot = piv_z.round(2).astype(str)
        for r in piv_z.index:
            for c in piv_z.columns:
                p = piv_p.loc[r, c]
                star = "**" if p < 0.01 else ("*" if p < 0.05 else ("+" if p < 0.10 else ""))
                annot.loc[r, c] = f"{piv_z.loc[r, c]:.2f}{star}"

        fig, ax = plt.subplots(figsize=(max(9, piv_z.shape[1] * 1.4), max(5, piv_z.shape[0] * 0.42)))
        sns.heatmap(
            piv_z,
            cmap="RdBu_r",
            center=0,
            vmin=-3.0,
            vmax=4.0,
            annot=annot,
            fmt="",
            linewidths=0.5,
            cbar_kws={"label": "scDRS Association MC Z-score (assoc_mcz; ** p<0.01, * p<0.05, + p<0.10)"},
            ax=ax,
        )
        ax.set_title(f"scDRS Disease Relevance Across {atlas_name} ({primary_gcol})", fontsize=12, fontweight="bold")
        ax.set_xlabel("Unbiased PSD Genetic Gene Set")
        ax.set_ylabel("Cell Type / Subcluster")
        plt.xticks(rotation=35, ha="right")
        plt.tight_layout()
        fig_path = FIG_DIR / f"scdrs_heatmap_{atlas_name}.png"
        plt.savefig(fig_path, dpi=300)
        plt.close(fig)
        print(f"Saved figure -> {fig_path}")

    # Plot 2: Gene specificity heatmap in Tabula Sapiens Skin and 12-Dataset Skin Atlas
    for spec_f in sorted(RESULTS_DIR.glob("gene_specificity_*.tsv")):
        atlas_name = spec_f.stem.replace("gene_specificity_", "")
        df_s = pd.read_csv(spec_f, sep="\t")
        if df_s.empty:
            continue
        idx_col = "cell_type" if "cell_type" in df_s.columns else "compartment"
        piv_s = df_s.pivot(index=idx_col, columns="gene", values="specificity_z").fillna(0.0)
        fig, ax = plt.subplots(figsize=(max(10, piv_s.shape[1] * 0.55), max(5, piv_s.shape[0] * 0.38)))
        sns.heatmap(
            piv_s,
            cmap="vlag",
            center=0,
            vmin=-2.5,
            vmax=2.5,
            linewidths=0.4,
            cbar_kws={"label": "Cell-Type Expression Specificity (Z-score)"},
            ax=ax,
        )
        ax.set_title(f"PSD Risk Gene Expression Specificity Across {atlas_name}", fontsize=12, fontweight="bold")
        ax.set_xlabel("PSD Candidate Risk Genes (GWAS + OpenTargets + UKB Exome)")
        ax.set_ylabel("Cell Type / Subcluster")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        fig_path = FIG_DIR / f"gene_specificity_heatmap_{atlas_name}.png"
        plt.savefig(fig_path, dpi=300)
        plt.close(fig)
        print(f"Saved figure -> {fig_path}")


if __name__ == "__main__":
    main()
